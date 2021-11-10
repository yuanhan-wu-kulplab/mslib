#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input residues
		required.push_back("pdb");
		required.push_back("nanoparticle");
		required.push_back("radius");
		required.push_back("positionNano");
		required.push_back("positionImmunogen");
		required.push_back("geocenterNanoPdb");
		required.push_back("rotate");
		required.push_back("linker");
		required.push_back("nanoFoldSymmetry");
		// Max C-alpha clashes allowed between potential fusions
		optional.push_back("maxCaClashes");



	}

	// Storage for the vales of each optional
	string pdb;
  	string nanoparticle;
        double maxCaClashes;
        double radius;
        string positionNano;
        string geocenterNano;
        vector<string> positionImmunogen;
        string linker;
        bool rotate;
        int foldSymmetry;
  
	// Storage for different types of options
	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);
Frame computeFrameFromTermini(System &_sys);
