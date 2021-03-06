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
    if(argc < 3)
    {
        cerr << "Need to pass:" << endl;
        cerr << "1) Path to folder containing the data" << endl;
        cerr << "2) g2 values file" << endl;
        return 1;
    }

    // Some declarations for later
    string prefix = "GEOM";
    string path = argv[1];
    string name = "ij2";



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
    
    // Read g2 values from file path/argv[2]
    string g2_filename = path + argv[2];

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
    

    // Open output file 
    string out_filename_A = path + "/observables/A" + name + ".txt";
    string out_filename_B = path + "/observables/B" + name + ".txt";
    ofstream out_obs_A;
    ofstream out_obs_B;
    out_obs_A.open(out_filename_A);
    out_obs_B.open(out_filename_B);

    if(!out_obs_A || !out_obs_B)
    {
        cerr << "Error: file " + out_filename_A + " or " + out_filename_B + " could not be opened." << endl;
        return 1;
    }

    // Cycle on g2 values
    for(const auto& g2 : g2_vec)
    {
        // Print value of g2 being processed
        clog << "g2: " << g2 << endl;

        // Create vector of uncorrelated samples
        vec samples_A(job_vec.size());
        vec samples_B(job_vec.size());

        // Cycle on jobs in the array
        for(unsigned i=0; i<job_vec.size(); ++i)
        {
            // Open data files
            string array_path = path + "/" + cc_to_name(g2) + "/" + to_string(job_vec[i]);
            string filename = data_to_name(sm.p, sm.q, sm.dim, g2, prefix);
            ifstream in_hl;
            in_hl.open(array_path + "/" + filename + "_HL.txt");

            // Create vector of correlated samples
            vec vec_corr_A(sm.samples);
            vec vec_corr_B(sm.samples);

            // Cycle on samples
            for(int j=0; j<sm.samples; ++j) 
            {
                Geom24 G(sm.p, sm.q, sm.dim, g2);
                G.read_mat(in_hl);

                cx_mat W = G.get_mat(0) + cx_double(0.,1.)*G.get_mat(1);
                cx_double trW = trace(W);
                
                cx_mat V = (abs(trW)/trW)*( W - (trW/double(sm.dim))*cx_mat(sm.dim, sm.dim, fill::eye) );
                cx_mat A = 0.5*(V + V.t());
                cx_mat B = 0.5*cx_double(0.,-1.)*(V - V.t());

                if(abs(trace(V)) > 1e-8)
                {
                    cerr << "Error: V is not traceless." << endl;
                    return 1;
                }
                if(!A.is_hermitian())
                {
                    cerr << "Error: A is not hermitian." << endl;
                    return 1;
                }
                if(!B.is_hermitian())
                {
                    cerr << "Error: B is not hermitian." << endl;
                    return 1;
                }


                // ***** COMPUTE OBSERVABLE HERE *****
                double temp_A = 0;
                double temp_B = 0;
                int counter = 0;
                for(int i=0; i<sm.dim; ++i)
                {
                    for(int j=i+1; j<sm.dim; ++j)
                    {
                        temp_A += (A(i,j)*A(j,i)).real();
                        temp_B += (B(i,j)*B(j,i)).real();

                        ++counter;
                    }
                }
                temp_A /= counter;
                temp_B /= counter;
                // ***** THAT'S IT, YOU'RE DONE *****

                vec_corr_A(j) = temp_A;
                vec_corr_B(j) = temp_B;
            }
            in_hl.close();

            // Initialize i-th element of vector of uncorrelated samples with mean of job #i
            samples_A(i) = mean(vec_corr_A);
            samples_B(i) = mean(vec_corr_B);
        }


        // Output mean and error of observable
        double avg = 0;
        double var = 0;
        jackknife(samples_A, avg, var, my_mean);
        double err = sqrt(var);
        out_obs_A << g2 << " " << avg << " " << err << endl;
        jackknife(samples_B, avg, var, my_mean);
        err = sqrt(var);
        out_obs_B << g2 << " " << avg << " " << err << endl;
    }

    out_obs_A.close();
    out_obs_B.close();

    //********* END ANALYSIS **********//

    return 0;
}
