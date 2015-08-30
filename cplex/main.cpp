#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include "cpxmacro.h"

using namespace std;


// CPlex globals
int status;
char errmsg[BUF_SIZE];

// Given the variables array with N*N y and N*N x
// Outputs the index of a given x_xy or y_xy
int idx(char v, int x, int y, int N) {
    if (v == 'x') {
        return x*N+y+N*N;
    }
    else {
        return x*N+y;
    }
}

int main (int argc, char const *argv[])
{
    try
    {
        
        /*
         * 
         * INIT VARIABILI NECESSARIE
         * 
         */
        
        DECL_ENV( env );
        DECL_PROB( env, lp );
        
        //TODO pls check argv[1] exists or die
        string problem_file = argv[1];
        ifstream data(problem_file.c_str());
        
        string line;
        vector<double>* weights = new vector<double>();
        vector<char>* variables = new vector<char>();
        vector<double>* lower_bounds = new vector<double>();
        vector<double>* upper_bounds = new vector<double>();
        int N = 0;
        
        
        /*
         * 
         * FUNZIONE OBIETTIVO
         * 
         */
        
        // Leggi i pesi, crea le variabili y del modello dato
        while(getline(data,line))
        {
            stringstream lineStream(line);
            string cell;
            while(getline(lineStream,cell,','))
            {
                double weight = atof(cell.c_str());
                weights->push_back(weight);
                variables->push_back('B');
                lower_bounds->push_back(0);
                upper_bounds->push_back(1);
            }
            N++;
        }
        
        int acard = weights->size();
        
        // Crea le variabili x del modello dato
        for (int i = 0; i < acard; i++) {
            weights->push_back(0);
            variables->push_back('I');
            lower_bounds->push_back(0);
            upper_bounds->push_back(N);        
        }

        char ** xname = NULL;

        // Funzione obiettivo, vincoli di tipo delle variabili
        CHECKED_CPX_CALL(CPXnewcols, env, lp, weights->size(), &(*weights)[0], &(*lower_bounds)[0], &(*upper_bounds)[0], &(*variables)[0], xname);
        
        
        
        /*
         * 
         * CREA I VINCOLI DEL PROBLEMA
         * 
         */

        
        // Inizializzo strutture dati per aggiungere constraints
        // numeric rhs
        double rhs[1] = {0};        
        // constraint type
        char sense[1] = {'E'};        
        // constraint coefficients
        vector<double>* rmatval = new vector<double>();
        // How many coeffs per row
        int rmatbeg[1] = {0};
        // Position of the constraint coefficient
        vector<int>* rmatind = new vector<int>();
        // Stuff
        char ** newcolnames = NULL;
        char ** rownames = NULL;
        
        
        // Vincolo 1
        
        for (int i=1; i<N; i++) {
            rmatval->push_back(1.0);
            rmatind->push_back(idx('x',0,i,N));
        }
        rhs[0] = N;
        CHECKED_CPX_CALL( CPXaddrows, env, lp, 0, 1, N-1, &rhs[0], &sense[0], &rmatbeg[0], &(*rmatind)[0], &(*rmatval)[0], newcolnames , rownames );
        
        rmatval->clear();
        rmatind->clear();
        
        
        // Vincolo 2
        
        for (int k=1; k<N; k++) {
            for (int i=0; i<N; i++) {
                rmatval->push_back(1);
                rmatind->push_back(idx('x',i,k,N));
                rmatval->push_back(-1);
                rmatind->push_back(idx('x',k,i,N));
            }
            rhs[0] = 1;
            CHECKED_CPX_CALL( CPXaddrows, env, lp, 0, 1, rmatval->size(), &rhs[0], &sense[0], &rmatbeg[0], &(*rmatind)[0], &(*rmatval)[0], newcolnames , rownames );
            rmatval->clear();
            rmatind->clear();
        }
        
        
        // Vincolo 3
        
        for (int i=0; i<N; i++) {
            for (int j=0; j<N; j++) {
                rmatval->push_back(1);
                rmatind->push_back(idx('y',i,j,N));
            }
            rhs[0] = 1;
            CHECKED_CPX_CALL( CPXaddrows, env, lp, 0, 1, rmatval->size(), &rhs[0], &sense[0], &rmatbeg[0], &(*rmatind)[0], &(*rmatval)[0], newcolnames , rownames );
            rmatval->clear();
            rmatind->clear();
        }
        
        
        // Vincolo 4
        
        for (int j=0; j<N; j++) {
            for (int i=0; i<N; i++) {
                rmatval->push_back(1);
                rmatind->push_back(idx('y',i,j,N));
            }
            rhs[0] = 1;
            CHECKED_CPX_CALL( CPXaddrows, env, lp, 0, 1, rmatval->size(), &rhs[0], &sense[0], &rmatbeg[0], &(*rmatind)[0], &(*rmatval)[0], newcolnames , rownames );
            rmatval->clear();
            rmatind->clear();
        }
        
        
        // Vincolo 5
        
        for (int i=0; i<N; i++) {
            for (int j=0; j<N; j++) {
                rmatval->push_back(1);
                rmatind->push_back(idx('x',i,j,N));
                rmatval->push_back(-N);
                rmatind->push_back(idx('y',i,j,N));
                
                rhs[0] = 0;
                sense[0] = 'L';
                CHECKED_CPX_CALL( CPXaddrows, env, lp, 0, 1, rmatval->size(), &rhs[0], &sense[0], &rmatbeg[0], &(*rmatind)[0], &(*rmatval)[0], newcolnames , rownames );
                rmatval->clear();
                rmatind->clear();
            }
        }
        
        
        // Ora rimuovo le colonne relative ai cappi
        int delstat[N*N*2];
        for (int i=0; i<N*N*2; i++)
            delstat[i] = 0;
        for (int i=0; i<N; i++) {
            delstat[idx('x',i,i,N)] = 1;
            delstat[idx('y',i,i,N)] = 1;
        }
        
        CHECKED_CPX_CALL( CPXdelsetcols, env, lp, &delstat[0] );


        /*
         * 
         * OTTIMIZZA
         * 
         */
        double start, end, total_time;
        CHECKED_CPX_CALL(CPXgettime, env, &start);
        CHECKED_CPX_CALL( CPXmipopt, env, lp );
        CHECKED_CPX_CALL(CPXgettime, env, &end);
        total_time = end - start;
        
        /*
         * 
         * STAMPA RISULTATI
         * 
         */
        
        double objval;
        CHECKED_CPX_CALL( CPXgetobjval, env, lp, &objval );
        string results_file = "../results/cplex/";
        string filename(argv[1]);
        results_file =  results_file + filename.substr(filename.rfind('/')+1);
        ofstream res (results_file.c_str());
        
        if (res.is_open())
        {
            res << fixed << setprecision(8) << objval << endl;
            res << fixed << setprecision(8) << total_time;
            res.close();
        }

        /*
         * 
         * CHIUDI
         * 
         */
        
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
    }
    catch(std::exception& e)
    {
        std::cout << ">>>EXCEPTION: " << e.what() << std::endl;
    }
    return 0;
}