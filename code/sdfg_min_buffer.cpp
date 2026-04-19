/*
Anusheel Nand
Find SAS that minimizes SDF's total buffer size
Dynamic Programming
EEC 284, Spring 2026
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <unordered_map>
#include <numeric>
#include <vector>
#include <limits>
#include <tuple>
#include <sstream>

// Easy tracking of rates per actor
struct rates {
    int produce;
    int consume;
};

// Track number of actors in graph for matrix sizing
int num_actors = 1;

char idx_to_actor (int idx);
std::map<char, rates> parse_sdfg(std::string fileName);
std::unordered_map<char, int> parse_sas(std::string fileName);
std::vector<std::vector<int>> fill_gcf (std::unordered_map<char, int> sas);
std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>>  min_buffer_sizes (std::map<char, rates> sdfg, 
                                                std::unordered_map<char, int> sas, const std::vector<std::vector<int>>& gcf);
std::string cuts_to_SAS(int L, int R, const std::vector<std::vector<int>>& gcf, const std::vector<std::vector<int>>& cuts);
void print_matrix (const std::vector<std::vector<int>>& matrix);

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cout << "Usage: ./min_buffer [sdfg file] [flat sas file]" << std::endl;
        exit(1);
    }

    std::map<char, rates> sdfg = parse_sdfg(argv[1]);
    std::unordered_map<char, int> flat_sas = parse_sas(argv[2]);

    // Print data flow graph from arg to terminal
    std::cout << "SDFG Graph:" << std::endl;
    int act_count = 0;
    for (auto actor : sdfg)
    {
        act_count++;
        if(actor.second.consume != 0)
        {
            std::cout << actor.second.consume << " ";
        }
        std::cout << "("  << actor.first << ")";
        if(act_count != num_actors)
        {
            std::cout << " " << actor.second.produce << "-->";
        }
        else
        {
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;

    // Print flat SAS from arg to terminal
    std::cout << "Flat SAS:" << std::endl;
    for(int i = 0; i < num_actors; i++)
    {
        char actor = idx_to_actor(i);
        std::cout << flat_sas[actor] << actor;
        if(i == num_actors-1)
        {
            std::cout << std::endl;
        }
        else
        {
            std::cout << " ";
        }
    }
    std::cout << std::endl;

    std::vector<std::vector<int>> gcf = fill_gcf(flat_sas);
    std::cout << "GCF Matrix: " << std::endl; 
    print_matrix(gcf);

    std::vector<std::vector<int>> sizes(num_actors, std::vector<int>(num_actors, 0));
    std::vector<std::vector<int>> cuts(num_actors, std::vector<int>(num_actors, 0));
    std::tie(sizes, cuts) = min_buffer_sizes(sdfg, flat_sas, gcf);

    std::cout << "Buffer Sizes Matrix: " << std::endl; 
    print_matrix(sizes);
    std::cout << "Cuts Matrix: " << std::endl; 
    print_matrix(cuts);

    std::string min_SAS = cuts_to_SAS(0, num_actors-1, gcf, cuts);

    std::cout << "SAS for minimum buffer: " << min_SAS << std::endl;
    std::cout << "Minimum buffer size: " << sizes[0][num_actors-1] << std::endl;

    return 0;
}

// Hacky way to convert from gcf/sizes matrix index to actor name
char idx_to_actor (int idx)
{
    return idx+65;
}

// Load the synchronous data flow graph (SDFG) from a file into a map
std::map<char, rates> parse_sdfg(std::string fileName)
{
    std::ifstream txt(fileName);
    if (txt.fail())
    {
        std::cout << "Invalid file" << std::endl;
        exit(1);
    }

    char start, end;
    int produce, consume;

    std::map<char, rates> sdfg;

    while (txt >> start >> produce >> consume >> end)
    {
        sdfg[start].produce = produce;
        sdfg[end].consume = consume;
        num_actors++;
    }
    txt.close();
    return sdfg;
}

// Load the flat single actor schedule from a file into hash table
std::unordered_map<char, int> parse_sas(std::string fileName)
{
    std::ifstream txt(fileName);
    if (txt.fail())
    {
        std::cout << "Invalid file" << std::endl;
        exit(1);
    }

    char actor;
    int occurences;

    std::unordered_map<char, int> flat_sas;

    while (txt >> actor >> occurences)
    {
        flat_sas[actor] = occurences;
    }
    txt.close();
    return flat_sas;
}

// Fill matrix of greatest common factors for easy access during DP loop
std::vector<std::vector<int>> fill_gcf (std::unordered_map<char, int> sas)
{
    std::vector<std::vector<int>> gcf(num_actors, std::vector<int>(num_actors, 0));
    for (int i = 0; i < num_actors; i++)
    {
        gcf[i][i] = sas[idx_to_actor(i)];
        for (int j = i+1; j < num_actors; j++)
        {
            int j_occurences = sas[idx_to_actor(j)];
            gcf[i][j] = std::gcd(gcf[i][j-1], j_occurences);
        }
    }
    return gcf;
}

// Dynamic programming buffer minimization

// Returns matrix containing minimal buffer sizes for range of actors
// Returns another matrix containing ideal cuts for later reconstruction of minimal buffer SAS
std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>>  min_buffer_sizes (std::map<char, rates> sdfg, std::unordered_map<char, int> sas, const std::vector<std::vector<int>>& gcf)
{
    std::vector<std::vector<int>> sizes(num_actors, std::vector<int>(num_actors, 0));
    std::vector<std::vector<int>> cuts(num_actors, std::vector<int>(num_actors, 0));
    // No need for single actor fill since sizes[][] matrix is initialized with zeroes
    for (int k = 1; k < num_actors; k++)
    {
        for(int i = 0; i < (num_actors - k); i++)
        {
            int j = i + k;
            int min_buffer_size = std::numeric_limits<int>::max();
            for(int cut = i; cut < j; cut++)
            {
                int buffer_size = sizes[i][cut] + sizes[cut+1][j] + 
                        ((sas[idx_to_actor(cut)] * sdfg[idx_to_actor(cut)].produce) / (gcf[i][j]));
                if(buffer_size < min_buffer_size)
                {
                    min_buffer_size = buffer_size;
                    sizes[i][j] = buffer_size;
                    cuts[i][j] = cut;
                }
            }
        }
    }

    return {sizes, cuts};
}

// Heavily based on reference's Figure 4: 
// Praveen K. Murthy, Shuvra S. Bhattacharyya, Edward A. Lee: Joint Minimization of Code and Data for Synchronous Dataflow Programs. 
// Formal Methods in System Design 11(1): 41-70 (1997)
// String formatting for the SAS was guided by Claude Opus 4.6
std::string cuts_to_SAS(int L, int R, const std::vector<std::vector<int>>& gcf, const std::vector<std::vector<int>>& cuts)
{
    if (L == R)
    {
        std::string actor(1,idx_to_actor(L));
        return actor;
    }

    // Best cut to factor from
    int cut = cuts[L][R];

    int i_l = (gcf[L][cut]) / (gcf[L][R]);
    int i_r = (gcf[cut+1][R]) / (gcf[L][R]);

    std::string left = cuts_to_SAS(L, cut, gcf, cuts);
    std::string right = cuts_to_SAS(cut+1, R, gcf, cuts);

    std::stringstream left_half;
    if (i_l > 1)
    {
        left_half << "(" << i_l << " " << left << ")";
    }
    else
    {
        left_half << left;
    }

    std::stringstream right_half;
    if (i_r > 1)
    {
        right_half << "(" << i_r << " " << right << ")";
    }
    else
    {
        right_half << right;
    }

    std::stringstream comb;
    comb << left_half.str() << " " << right_half.str();
    return comb.str();
}

// Helper to print progress on matrices
void print_matrix (const std::vector<std::vector<int>>& matrix)
{
    std::cout << "  "; // Space to label left side of matrix
    for(int i = 0; i < matrix.size(); i++)
    {
        std::cout << std::setw(2) << idx_to_actor(i) << " ";
    }
    std::cout << std::endl;

    for(int i = 0; i < matrix.size(); i++)
    {
        std::cout << idx_to_actor(i) << " ";
        for (int j = 0; j < matrix[i].size(); j++)
        {
            std::cout << std::setw(2) << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}