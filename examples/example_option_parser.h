
using namespace std;

// DECLARE an Options structure (this usually goes in a header file like 'example_option_parser.h')
struct Options {

	// Set up options here...
	Options(){
		// PDB list
		required.push_back("pdb");
	}

	string pdb;

	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};

// DECLARE a setupOptions function to keep the main tidy (this usually goes in a header file like 'example_option_parser.h')
Options setupOptions(int theArgc, char * theArgv[]);
