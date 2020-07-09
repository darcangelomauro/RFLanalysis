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
        cerr << "1) Name of the folder containing the data" << endl;
        cerr << "2) First index of the jobs array" << endl;
        cerr << "3) Number of jobs in the array" << endl;
        cerr << "4) Name of the observable" << endl;
        return 1;
    }

    // Some declarations for later
    string prefix = "GEOM";
    string path = argv[1];
    int fst_jarr = stoi(argv[2]);
    int num_jarr = stoi(argv[3]);
    string name = argv[4]; 



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

    if(sm.p!=1 || sm.q!=3)
    {
        cerr << "Error: geometry is not (1,3)" << endl;
        return 1;
    }
    //********* END PARAMETER INITIALIZATION **********//



    //********* BEGIN ANALYSIS **********//
    

    // Open output file 
    string out_filename = path + "/observables/" + name + ".txt";
    ofstream out_obs;
    out_obs.open(out_filename);

    if(!out_obs)
    {
        cerr << "Error: file " + out_filename + " could not be opened." << endl;
        return 1;
    }

    // Cycle on g2 values
    double g2 = sm.g2_i;
    while(g2 < sm.g2_f)
    {
        // Print value of g2 being precessed
        clog << "g2: " << g2 << endl;

        // Create vector of uncorrelated samples
        vec samples(num_jarr);

        // Cycle on jobs in the array
        for(int i=0; i<num_jarr; ++i)
        {
            // Open data files
            string array_path = path + "/" + to_string(i+fst_jarr);
            string filename = filename_from_data(sm.p, sm.q, sm.dim, g2, prefix);
            ifstream in_s, in_hl;
            in_s.open(array_path + "/" + filename + "_S.txt");
            in_hl.open(array_path + "/" + filename + "_HL.txt");

            // Create vector of correlated samples
            int length = n_meas(sm.iter_simul, sm.gap);
            vec vec_corr(length);

            // Cycle on samples
            for(int j=0; j<length; ++j) 
            {
                Geom24 G(sm.p, sm.q, sm.dim, g2);
                double S2, S4;
                in_s >> S2 >> S4;
                G.read_mat(in_hl);

                // ***** COMPUTE OBSERVABLE HERE *****
                cx_mat C = G.get_mat(2)*G.get_mat(5) - G.get_mat(5)*G.get_mat(2);
                cx_mat Ct = C.t();
                double temp = trace(Ct*C).real();
                // ***** THAT'S IT, YOU'RE DONE *****

                vec_corr(j) = temp;
            }
            in_s.close();
            in_hl.close();

            // Initialize i-th element of vector of uncorrelated samples with mean of job #i
            samples(i) = mean(vec_corr);
        }


        // Output mean and error of observable
        double avg = 0;
        double var = 0;
        jackknife(samples, avg, var, my_mean);
        double err = sqrt(var);
        out_obs << g2 << " " << avg << " " << err << endl;

        g2 += sm.g2_step;
    }

    out_obs.close();

    //********* END ANALYSIS **********//

    return 0;
}
