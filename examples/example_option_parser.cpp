
#include "OptionParser.h"

#include <iostream>
#include <vector>
#include <string>

#include "example_option_parser.h"

using namespace std;
using namespace MSL;

// Very simple main 
int main(int argc, char *argv[]) {

    // Do all option parsing with argc and argv
    Options opt = setupOptions(argc, argv);

    // Now use the options:
    cout << "\n\n***EXAMPLE START***\n\n"<< "PDB FILE IS: "<<opt.pdb<<"\n\n***EXAMPLE COMPLETE***\n\n";
}


// Must define this function
Options setupOptions(int theArgc, char * theArgv[]){

    // Create the options
    Options opt;

    // Parse the options
    OptionParser OP;
    OP.setRequired(opt.required);	
    OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option
    OP.readArgv(theArgc, theArgv);

    if (OP.countOptions() == 0){
	cout << "Usage: example_option_parser " << endl;
	cout << endl;
	cout << "\n";
	cout << "pdb PDBFILE\n";
	cout << endl;
	exit(0);
    }

    // Parse the 'pdb' option
    opt.pdb = OP.getString("pdb");
    if (OP.fail()){
	cerr << "ERROR 1111 no pdb specified."<<endl;
	exit(1111);
    }

    return opt;

}

