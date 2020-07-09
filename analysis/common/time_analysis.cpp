#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <armadillo>
#include "geometry.hpp"
#include "utils.hpp"
#include "params.hpp"
#include "statistics.hpp"

using namespace std;
using namespace arma;

int main(int argc, char** argv)
{
    // Check arguments
    if(argc < 5)
    {
        cerr << "Need to pass:" << endl;
        cerr << "1) Name of the folder containing the data" << endl;
        cerr << "2) First index of the jobs array" << endl;
        cerr << "3) Number of jobs in the array" << endl;
        cerr << "4) Slurm ID" << endl;
        return 1;
    }

    // Some declarations for later
    string path = argv[1];
    int fst_jarr = stoi(argv[2]);
    int num_jarr = stoi(argv[3]);
    string job_id = argv[4];


    //********* BEGIN ANALYSIS **********//
    

    // Open output file 
    string out_filename = path + "/observables/time.txt";
    ofstream out_obs;
    out_obs.open(out_filename);

    if(!out_obs)
    {
        cerr << "Error: file " + out_filename + " could not be opened." << endl;
        return 1;
    }

    // Create vector of uncorrelated samples
    Col<unsigned long> samples(num_jarr);

    // Cycle on jobs in the array
    for(int i=0; i<num_jarr; ++i)
    {
        // Read stdout from file path/stdout/slurm-(job_id)_(i+fst_jarr).out
        string in_filename = path + "/stdout/slurm-" + job_id + "_" + to_string(i+fst_jarr) + ".out";
        ifstream in;
        in.open(in_filename);

        if(in)
        {
            // Read until you reach the first timestamp
            string temp;
            vector<unsigned long> vec;
            while(in >> temp)
            {
                if(temp == "timestamp:")
                {
                    in >> temp;
                    vec.push_back(stoul(temp));
                }
            }

            
            unsigned long start = *(vec.begin());
            unsigned long end = *(vec.end()-1);
            
            // Initialize i-th element of vector of samples with time of job #i
            samples(i) = end-start;

            out_obs << i << " " << end-start << endl;
        
            in.close();
        }
        else
        {
            cout << "Error in opening stdout file: " + in_filename << endl;
            return 1;
        }
    }

    // Output mean and error of observable
    unsigned long avg = mean(samples);
    unsigned long err = stddev(samples)/sqrt(num_jarr);
    
    out_obs << "Average (seconds): " << avg << " " << err << endl;
    out_obs.close();

    cout << "Average time: " << avg/3600. << " +- " << err/3600. << " hours" << endl;

    //********* END ANALYSIS **********//

    return 0;
}
