#include <iostream>
#include <string>
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
        cerr << "1) Path to folder containing the data" << endl;
        cerr << "2) Name of the observable" << endl;
        cerr << "3) Coupling constant value" << endl;
        cerr << "4) Positive extremum of histogram" << endl;
        return 1;
    }

    // Some declarations for later
    string prefix = "GEOM";
    string path = argv[1];
    string name = argv[2]; 
    double g2_input = stod(argv[3]);
    double extr = abs(stod(argv[4]));



    //********* BEGIN PARAMETER INITIALIZATION **********//
    
    // Read simulation parameters from file path/init.txt
    string init_filename = path + "/init.txt";

    struct Simul_params sm;
    ifstream in_init;
    in_init.open(init_filename);

    if(!read_init_stream(in_init, sm))
    {
        cerr << "Error: couldn't read file " + init_filename << endl;
        return 1;
    }

    cout << "File " + init_filename + " contains the following parameters:" << endl;
    cout << sm.control << endl;

    if(!params_validity(sm))
    {
        cerr << "Error: file " + init_filename + " does not contain the necessary parameters." << endl;
        return 1;
    }

    in_init.close();

    //********* END PARAMETER INITIALIZATION **********//

    
    //********* BEGIN G2 INITIALIZATION **********//
    
    // Read simulation parameters from file path/init.txt
    string g2_filename = path + "/g2_val.txt";

    ifstream in_g2;
    in_g2.open(g2_filename);

    if(!in_g2.is_open())
    {
        cerr << "Error: couldn't read file " + g2_filename << endl;
        return 1;
    }

    vector<double> g2_vec;
    double temp_g2;
    while(in_g2 >> temp_g2)
        g2_vec.push_back(temp_g2);
    
    cout << "File " + g2_filename + " contains " << g2_vec.size() << " g2 values:" << endl;
    cout << "From " << *g2_vec.begin() << " to " << *(g2_vec.end()-1) << endl;

    in_g2.close();

    //********* END G2 INITIALIZATION **********//

    
    //********* BEGIN JOB ARRAY INITIALIZATION **********//
    
    // Read simulation parameters from file path/init.txt
    string job_filename = path + "/job_idx.txt";

    ifstream in_job;
    in_job.open(job_filename);

    if(!in_job.is_open())
    {
        cerr << "Error: couldn't read file " + job_filename << endl;
        return 1;
    }

    vector<int> job_vec;
    int temp_job;
    while(in_job >> temp_job)
        job_vec.push_back(temp_job);
    
    cout << "File " + job_filename + " contains " << job_vec.size() << " job indices:" << endl;
    cout << "From " << *job_vec.begin() << " to " << *(job_vec.end()-1) << endl;

    in_job.close();

    //********* END JOB ARRAY INITIALIZATION **********//


    
    //********* BEGIN ANALYSIS **********//
    


    // Cycle on g2 values
    for(const auto& g2 : g2_vec)
    {
        if(abs(g2-g2_input) < 1e-8)
        {
            // Print value of g2 being processed
            clog << "g2: " << g2 << endl;
                
            // Open output file 
            string out_filename = path + "/observables/" + name + "_" + cc_to_name(g2) + ".txt";
            ofstream out_obs(out_filename);
                
            if(!out_obs)
            {
                cerr << "Error: file " + out_filename + " could not be opened." << endl;
                return 1;
            }

            // Create vector that stores every eigenvalue of D of every sample of every job
            Geom24 T(sm.p, sm.q, 1, 1);
            int evals_num = sm.dim*sm.dim*T.get_dim_omega()*sm.samples*job_vec.size();
            vec evals(evals_num);
            clog << "Total number of eigenvalues for histogram: " << evals_num << endl;

            // Create a counter that tells me how many eigenvalues have been stored so far
            unsigned counter = 0;


            // Cycle on jobs in the array
            for(unsigned i=0; i<job_vec.size(); ++i)
            {
                clog << "job: " << job_vec[i] << endl;

            
                // Open data files
                string array_path = path + "/" + cc_to_name(g2) + "/" + to_string(job_vec[i]);
                string filename = data_to_name(sm.p, sm.q, sm.dim, g2, prefix);
                ifstream in_hl;
                in_hl.open(array_path + "/" + filename + "_HL.txt");

                // Cycle on samples
                for(int j=0; j<sm.samples; ++j) 
                {
                    Geom24 G(sm.p, sm.q, sm.dim, g2);
                    G.read_mat(in_hl);

                    // ***** COMPUTE OBSERVABLE HERE *****
                    cx_mat D = G.build_dirac();

                    vec temp = eig_sym(D);
                    for(const auto& val : temp)
                    {
                        evals(counter) = val;
                        ++counter;
                    }
                    // ***** THAT'S IT, YOU'RE DONE *****

                }
                in_hl.close();
            }

            if(counter != evals.n_elem)
            {
                cerr << "Error: number of samples not correct" << endl;
                return 1;
            }

            vec bins = linspace<vec>(-extr,extr,100);
            uvec dos = hist(evals, bins);

            for(unsigned i=0; i<dos.n_elem; ++i)
                out_obs << bins(i) << " " << double(dos(i))/evals_num << endl;

            out_obs.close();
        }
    }

    //********* END ANALYSIS **********//

    return 0;
}
