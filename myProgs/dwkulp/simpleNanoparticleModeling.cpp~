#include <iostream>
#include <cstdlib>
#include <queue>
#include <algorithm>
#include <map>

#include "System.h"
#include "OptionParser.h"
#include "Frame.h"
#include "MslOut.h"
#include "SysEnv.h"
#include "SasaCalculator.h"
#include "PhiPsiStatistics.h"
#include "MonteCarloManager.h"
#include "Quench.h"
#include "PDBTopology.h"
#include "AtomSelection.h"
#include "PyMolVisualization.h"
#include "AtomSelection.h"
#include "AtomContainer.h"
#include "RandomNumberGenerator.h"
#include "PDBWriter.h"
#include "Helanal.h"

#include "simpleNanoparticleModeling.h"

using namespace std;
using namespace MSL;


// MslOut 
static MslOut MSLOUT("simpleNanoparticleModeling");
static SysEnv SYSENV;

  
int main(int argc, char *argv[]) {

  // Parse commandline options
  Options opt = setupOptions(argc,argv);

  // Read in original pdb
  System sys;
  sys.readStructureFile(opt.pdb);
  sys.getAtomPointers().saveCoor("orig");

  Frame cterFrame;
  AtomPointerVector immunogenAts;
  if (opt.positionImmunogen.size() == 0){
 
    // Use average of C-term and Geomteric Center to build an axis to align to global Z-axis, put ave_c_term on 0,0,0
    cterFrame = computeFrameFromTermini(sys);
    
  } else {

    // Use specified positions to create a PCA frame
    for (uint p = 0;p < opt.positionImmunogen.size();p++){

      if (!sys.positionExists(opt.positionImmunogen[p])){
	cout << "Skipping position: "<<opt.positionImmunogen[p]<<" does not exist."<<endl;
	continue;
      }

      if (sys.getPosition(opt.positionImmunogen[p]).atomExists("CA")){
	immunogenAts.push_back(&sys.getPosition(opt.positionImmunogen[p]).getAtom("CA"));
      }
    }

    cterFrame.computeFrameFromPCA(immunogenAts);
  }
  
  cterFrame.transformToGlobalBasis(sys.getAtomPointers());
  sys.getAtomPointers().saveCoor("center");
  sys.writePdb("immunogen_center.pdb");


  // Read in original pdb
  System nano;
  nano.readStructureFile(opt.nanoparticle);
  
  // Collect 'nanoPoints' or junctional points for trimer
  CartesianPoint nanoGC(0,0,0);
  vector<CartesianPoint> nanoPoints;
  
  if (opt.geocenterNano != ""){
    System gc;
    gc.readStructureFile(opt.geocenterNano);
    nanoGC = gc.getAtom(0).getCoor();
  } else {
      cout << "NANO CHAIN SIZE: "<<nano.chainSize()<<endl;
      for (uint c = 0; c < nano.chainSize();c++){
    
	//Position &cpos = nano.getChain(c).getPosition(nano.getChain(c).positionSize()-1);
	if (!nano.getChain(c).positionExists(opt.positionNano)){
	  cout << "ERROR chain "<<nano.getChain(c).getChainId()<<" does not have position: "<<opt.positionNano<<endl;
	  exit(1231);
	}
    
	Position &cpos = nano.getChain(c).getPosition(opt.positionNano);
	cout << "SELECTING POSITION: "<<cpos<<endl;
	CartesianPoint thisChainGC(cpos.getAtom("CA").getCoor());
	nanoPoints.push_back(thisChainGC);
	nanoGC += thisChainGC;
      }
      nanoGC /= nano.chainSize();
  }

  


  Transforms tm;
  nano.writeCif(MslTools::stringf("nano_orig_%s.cif",MslTools::getFileName(opt.nanoparticle).c_str()));
  nano.getAtomPointers().saveCoor("orig");
  tm.translate(nano.getAtomPointers(),nanoGC*-1.0);
  nano.writeCif(MslTools::stringf("nano_center_%s.cif",MslTools::getFileName(opt.nanoparticle).c_str()));
  
  
  PyMolVisualization pymol;
  nanoPoints.clear();
  for (uint c = 0; c < nano.chainSize();c++){
    Position &cpos = nano.getChain(c).getPosition(opt.positionNano);
    CartesianPoint thisChainGC(cpos.getAtom("CA").getCoor());
    nanoPoints.push_back(thisChainGC);
  }

  vector<CartesianPoint> nanoPointsSymmetry;
  if (opt.foldSymmetry > 1){
    
    // Compute which chains at this position are closest
    map<int,bool> nanoPointsUsed;
    for (uint a = 0; a < nanoPoints.size();a++){

      // Don't use this nanoPoint if already been used.
      if (nanoPointsUsed.find(a) != nanoPointsUsed.end()) continue;

      // Mark this nanoPoint as used
      nanoPointsUsed[a] = true;
      
      priority_queue<pair<double,int> > closestChains;
      for (uint b = 0; b < nanoPoints.size();b++){
	if (b == a) continue;
	double dist = nanoPoints[a].distance(nanoPoints[b]);
	closestChains.push(pair<double,int>(-dist,b));
      }

      int c = 1;
      CartesianPoint thisChainGC(nanoPoints[a]);
      while (closestChains.size() > 0 && c < opt.foldSymmetry){
	cout << "A: "<<nano.getChain(a).getChainId()<<" "<<closestChains.top().first<<" "<<nano.getChain(closestChains.top().second).getChainId()<<endl;

	// Make a geometric center from the closest nanoPoints
	thisChainGC += nanoPoints[closestChains.top().second];

	// Mark these neighboring nano points as used
	nanoPointsUsed[closestChains.top().second] = true;
	
	closestChains.pop();
	c++;
      }
      thisChainGC = thisChainGC / opt.foldSymmetry;
      
      nanoPointsSymmetry.push_back(thisChainGC);
      
    }
    nanoPoints = nanoPointsSymmetry;
  } 

  cout << "NanoPoints: "<<nanoPoints.size()<<endl;
  
  // IF OPT.LINKER: Read in linker, compute frame
  System link;
  Frame linkerFrame;
  if (opt.linker != ""){
    link.readStructureFile(opt.linker);
    link.saveCoor("orig");
    linkerFrame.computeFrameFromPCA(link.getAtomPointers());
  }
  
  vector<string> names;
  vector<CartesianPoint> gcs;


  
  int totalClashes = 0;
  for (uint n = 0; n < nanoPoints.size();n++){
    cout << "Subunit: "<<n<<endl;
    pymol.createAtom(nanoPoints[n],MslTools::stringf("nanoPoint_%03d",n));
    
    // Compute a frame for trimer at surface 
    tm.translate(sys.getAtomPointers(),CartesianPoint(0,0,nanoPoints[n].distance(CartesianPoint(0,0,0))));
    tm.align(sys.getAtomPointers(),CartesianPoint(0,0,1),nanoPoints[n]);
    
    Frame surfaceFrame;
    if (opt.positionImmunogen.size() == 0) {
      surfaceFrame = computeFrameFromTermini(sys);
    } else {
      surfaceFrame.computeFrameFromPCA(immunogenAts);
    }
    
//    pymol.createFrame(surfaceFrame, MslTools::stringf("frameAtNanoSurface_%02d",n));
//    
//    ofstream fout;
//    fout.open(MslTools::stringf("trimerAtNanoSurface.py"));
//    fout << pymol;
//    fout.close();

    PDBWriter pout;
    pout.open(MslTools::stringf("trimerAtNanoSurface-%02d.pdb",n));
    pout.write(sys.getAtomPointers());
    pout.close();

    sys.saveCoor("currentSurface");
    sys.applySavedCoor("center");
    
    tm.translate(sys.getAtomPointers(),CartesianPoint(0,0,nanoPoints[n].distance(CartesianPoint(0,0,0))+opt.radius));
    tm.align(sys.getAtomPointers(),CartesianPoint(0,0,1),nanoPoints[n]);
    
    string name = MslTools::stringf("%s_trimer_n%03d_r%03d.pdb",MslTools::getFileName(opt.nanoparticle).c_str(),n,(int)opt.radius);
    sys.getAtomPointers().saveCoor(name);
    sys.writePdb(name);

    if (!opt.rotate) {
      sys.applySavedCoor("center");
      continue;
    }

    if (n == 0){
	string name2 = MslTools::stringf("%s_trimer_ROTATE_n%03d_r%03d.pdb",MslTools::getFileName(opt.nanoparticle).c_str(),n,(int)opt.radius);
	sys.getAtomPointers().saveCoor(name2);
	sys.writePdb(name2);
	names.push_back(name2);
	gcs.push_back(sys.getAtomPointers().getGeometricCenter());
	
	sys.applySavedCoor("center");
	continue;
    }




    
    // Create a copy of the current atoms
    AtomContainer currentAts(sys.getAtomPointers());
    currentAts.saveCoor("current");

    // Select the CA atoms from the copy
    AtomSelection sel(currentAts.getAtomPointers());
    AtomPointerVector ats = sel.select("name CA");	

    // Select the CA atoms from the system
    AtomSelection sel2(sys.getAtomPointers());
    AtomPointerVector ats2 = sel2.select("name CA");
	
    cout << "Rotating.."<<endl;

    /*
       Rotating algorithm:
         1. Random rotations in relative x,y,z (X -> upto 30 degrees, Y -> upto 30 degrees, Z -> upto 360 degrees)
	 2. Check for clashes
	 3. Save coordinates if low clashes
	 4. At end align sys "center" to lowest clashing conformation
     */
    
    RandomNumberGenerator rng;
    int minClashes = MslTools::intMax;
    int numRotSteps = 5000;
    string minString = "";
    for (uint rots = 0; rots < numRotSteps; rots++){
      
      // Select random angles in appropriate ranges..
      int angX = rng.getRandomInt(0,30)*rng.getRandomBit();
      int angY = rng.getRandomInt(0,30)*rng.getRandomBit();
      int angZ = rng.getRandomInt(120);

      // If 3/4 done, last quarter allow for large rotations in x and y
      if (rots > 0.75*numRotSteps){
	angX = rng.getRandomInt(30,75)*rng.getRandomBit();
	angY = rng.getRandomInt(30,75)*rng.getRandomBit();
      }

      // Might be good to check CHORD distance from previous position to see if this is reasonable. If you move 5x the bounding box it is probably not realistic.
      
      cout << "TRY: "<<rots<<" ROTATE: "<<angX<<","<<angY<<","<<angZ<<endl;
      // Rotate
      //PDBWriter pout1;
      //pout1.open("test_trimer1.pdb");
      //pout1.write(currentAts.getAtomPointers());
      //pout1.close();
      
      tm.rotate(currentAts.getAtomPointers(),angX,surfaceFrame["X"].getCenter()-surfaceFrame["X"].getDirection(),surfaceFrame["X"].getCenter());
      tm.rotate(currentAts.getAtomPointers(),angY,surfaceFrame["Y"].getCenter()-surfaceFrame["Y"].getDirection(),surfaceFrame["Y"].getCenter());
      tm.rotate(currentAts.getAtomPointers(),angZ,surfaceFrame["Z"].getCenter()-surfaceFrame["Z"].getDirection(),surfaceFrame["Z"].getCenter());

      //CartesianPoint xcen = surfaceFrame["X"].getCenter();
      //CartesianPoint xdir = surfaceFrame["X"].getDirection();
      //CartesianPoint ycen = surfaceFrame["Y"].getCenter();
      //CartesianPoint ydir = surfaceFrame["Y"].getDirection();
      //CartesianPoint zcen = surfaceFrame["Z"].getCenter();
      //CartesianPoint zdir = surfaceFrame["Z"].getDirection();
      //pymol.createArrow(xcen,xdir,"Xaxis");
      //pymol.createArrow(ycen,ydir,"Yaxis");
      //pymol.createArrow(zcen,zdir,"Zaxis");

      
      //ofstream fout;
      //fout.open(MslTools::stringf("test.py"));
      //fout << pymol;
      //fout.close();
      //
      //PDBWriter pout;
      //pout.open("test_trimer2.pdb");
      //pout.write(currentAts.getAtomPointers());
      //pout.close();
      //exit(0);
      // Get geometric center of ats
      CartesianPoint gc1 = ats.getGeometricCenter();
      
      // Do clash check
      int numClashes = 0;
      for (uint s = 0; s < names.size();s++){
	  
	//cout << "Clash check against name: "<<names[s]<<endl;

	CartesianPoint gc2 = gcs[s];
	double dist2 = gc1.distance2(gc2);

	// A neighbor is way to close < 55 Angstroms
	if (dist2 < 3025){
	  cout << "Number of clashes between "<<name<<" and "<<names[s]<< " is too close."<<endl;
	  numClashes = MslTools::intMax;
	  break;
	}

	// Test clashes for trimers less than 150 Angstroms
	if (dist2 < 22500){
	  
	  sys.applySavedCoor(names[s]);
	  
	  // Compute numClashes
	  int numPairClashes = 0;
	  for (uint a1 = 0; a1 < ats.size();a1++){
	    for (uint a2 = 0; a2 < ats2.size();a2++){
	      if (ats(a1).distance2(ats2(a2)) < 625){
		numPairClashes++;
	      }
	    }
	  }
	cout << "Number of clashes between "<<name<<" and "<<names[s]<< " is "<<numPairClashes<<" GC: "<<ats.getGeometricCenter().distance(ats2.getGeometricCenter())<<endl;

	numClashes += numPairClashes;
	}
	
      }
	
      if (numClashes < minClashes){
	currentAts.saveCoor("minClashes");
	minClashes = numClashes;
	cout << "New clashes: "<<minClashes<<" "<<angX<<" "<<angY<<" "<<angZ<<endl;
	minString = MslTools::stringf("MinClashes: %6d, Angles: %3d, %3d, %3d, Rotation: %3d",minClashes,angX,angY,angZ,rots);
      }

      if (numClashes == 0){
	minClashes = 0;
	minString = MslTools::stringf("MinClashes: %6d, Angles: %3d, %3d, %3d, Rotation: %3d",minClashes,angX,angY,angZ,rots);
	break;
      }

      // Reset the trimer
      currentAts.applySavedCoor("current");	
    }
    cout << "FINAL CLASH COUNT for subunit "<<n<< " is "<<minString<<endl;
    
    totalClashes += minClashes;
    
    // Align sys to currentAts
    tm.rmsdAlignment(ats2,ats,sys.getAtomPointers());
    cout <<"RMSD of "<<n<<" is: "<<tm.getLastRMSD()<<endl;
    // Store this sys orientation and write a pdb of it
    string name2 = MslTools::stringf("%s_trimer_ROTATE_n%03d_r%03d.pdb",MslTools::getFileName(opt.nanoparticle).c_str(),n,(int)opt.radius);
    sys.saveCoor(name2);
    sys.writePdb(name2);

    // Store name for next iteration
    names.push_back(name2);
    gcs.push_back(sys.getAtomPointers().getGeometricCenter());
    
    // Linker modeling algorithm:
    //    align line (frame using PCA from CA atoms?) to surface frame to current frame of trimer 
    if (opt.linker != ""){

      // Compute a frame from the current position of sys
      Frame currentFrame; 
      if (opt.positionImmunogen.size() == 0) {
	currentFrame = computeFrameFromTermini(sys);
      } else {
	currentFrame.computeFrameFromPCA(immunogenAts);
      }


	  
      // Compute a line between the nanoparticle surface and the trimer
      Line trimerLine((surfaceFrame["X"].getCenter()+currentFrame["X"].getCenter())/2,currentFrame["X"].getCenter()-surfaceFrame["X"].getCenter());

      trimerLine.setName("trimerLine");
      trimerLine.setOutputFormat("pymol");
      ofstream fout;
      fout.open(MslTools::stringf("test_%d.py",n));
      fout << trimerLine.toString();
      linkerFrame.setName("linkerFrame1");
      fout << linkerFrame.toString();


      // Translate frame to center of trimer line
      tm.translate(link.getAtomPointers(),trimerLine.getCenter()-linkerFrame["X"].getCenter());

      cout << "LINKER 2"<<endl;
      linkerFrame.computeFrameFromPCA(link.getAtomPointers());
      linkerFrame.setName("linkerFrame2");
      fout << linkerFrame.toString();
      
      // Compute an angle between Z-axis (long axis) to trimer Line
      double angle = linkerFrame["Z"].angle(trimerLine);

      // Compute the axis as a cross product so we can get lines to align
      CartesianPoint axis(linkerFrame["Z"].getDirection().cross(trimerLine.getDirection()));
      axis = axis.getUnit();
      tm.rotate(link.getAtomPointers(),-angle,linkerFrame["Z"].getCenter()-axis,linkerFrame["Z"].getCenter());


      linkerFrame.computeFrameFromPCA(link.getAtomPointers());
      linkerFrame.setName("linkerFrame3");
      fout << linkerFrame.toString();
      fout.close();


      // Write out linker PDB
      string nameLinker = MslTools::stringf("%s_linker_ROTATE_n%03d_r%03d.pdb",MslTools::getFileName(opt.nanoparticle).c_str(),n,(int)opt.radius);
      link.writePdb(nameLinker);

      // Reset linker coords + linkerFrame
      link.applySavedCoor("orig");
      linkerFrame.computeFrameFromPCA(link.getAtomPointers());
      
    }
    
    // Clean up the copied system
    currentAts.removeAllAtoms();

    // Reset system atoms to center
    sys.applySavedCoor("center");
  }
  cout << "Tally of clashes:  "<<totalClashes<<endl;

  


  
  // Compute total clashes of final model
  AtomSelection sel2(sys.getAtomPointers());
  AtomPointerVector ats2 = sel2.select("name CA");
  
  int numClashes = 0;
  for (uint s = 0; s < names.size();s++){

    sys.applySavedCoor(names[s]);
    AtomContainer currentAts(sys.getAtomPointers());
    AtomSelection sel(currentAts.getAtomPointers());
    AtomPointerVector ats = sel.select("name CA");
    
    for (uint s2 = s+1;s2 < names.size();s2++){
      sys.applySavedCoor(names[s2]);

    
      // Compute numClashes
      int numPairClashes = 0;
      for (uint a1 = 0; a1 < ats.size();a1++){
	for (uint a2 = 0; a2 < ats2.size();a2++){
	  if (ats(a1).distance2(ats2(a2)) < 625){
	    numPairClashes++;
	  }
	}
      }
      numClashes += numPairClashes;
      if (numPairClashes > 0){
	cout << "CLASHES BETWEEN: "<<names[s]<< " and "<<names[s2]<<" is "<<numPairClashes<<endl;
      }
      
    }

    currentAts.removeAllAtoms();
  }

  cout << "NUM TOTAL CLASHES: "<<numClashes<<" on-the-fly clash count: "<<totalClashes<<endl;

  
  /*
    
  // Spiral sampled sphere code...

  CartesianPoint center(0,0,0);
  PyMolVisualization pymol;
  pymol.createAtom(center, "Center");
  
  double radius = 10;
  double num    = 24;

  for (uint i = 0; i < num;i++){
    double phi = acos(1-2*i/num);
    double theta = M_PI * (1+sqrt(5) *i);
    CartesianPoint cp (radius*sin(phi)*cos(theta),radius*sin(theta)*sin(phi),radius*cos(phi));
    pymol.createAtom(cp, MslTools::stringf("CP-%i",i));
  }
  */
  
  ofstream fout;
  fout.open(MslTools::stringf("simpleNano.py"));
  fout << pymol;
  fout.close();
  
  MSLOUT.stream() << "Done."<<endl;
}


Options setupOptions(int theArgc, char * theArgv[]){
  Options opt;

  OptionParser OP;


  OP.setRequired(opt.required);
  OP.setAllowed(opt.optional);
  OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"
  OP.readArgv(theArgc, theArgv);

  if (OP.countOptions() == 0){
    cout << "Usage:" << endl;
    cout << endl;
    cout << "simpleNanoparticleModeling --pdb pdb\n";

    cout << "\nprogram options: "<<endl;
    for (uint i = 0; i < opt.required.size();i++){
      cout <<"R  --"<<opt.required[i]<<"  "<<endl;
    }
    cout <<endl;
    for (uint i = 0; i < opt.optional.size();i++){
      cout <<"O  --"<<opt.optional[i]<<"  "<<endl;
    }
    cout << endl;
    exit(0);
  }

  opt.pdb = OP.getString("pdb");
  if (OP.fail()){
    cerr << "ERROR 1111 pdb/cif not specified.\n";
    exit(1111);
  }

  opt.nanoparticle = OP.getString("nano");
  if (OP.fail()){
    cerr << "ERROR 1111 nano not specified.\n";
    exit(1111);
  }
  opt.radius = OP.getDouble("radius");
  if (OP.fail()){
    cerr << "ERROR 1111 radius not specified.\n";
    exit(1111);
  }

  opt.positionNano = OP.getString("positionNano");
  if (OP.fail()){
    cerr << "ERROR 1111 positionNano not specified.\n";
    exit(1111);
  }

  opt.positionImmunogen = OP.getStringVectorJoinAll("positionImmunogen");

  opt.geocenterNano = OP.getString("geocenterNanoPdb");
  opt.rotate = OP.getBool("rotate");
  opt.linker = OP.getString("linker");
  opt.foldSymmetry = OP.getInt("nanoFoldSymmetry");
  
  MSLOUT.stream() << "Options:\n"<<OP<<endl;
  return opt;
}





Frame computeFrameFromTermini(System &_sys){

  vector<CartesianPoint> points;
  
  for (uint c = 0; c < _sys.chainSize();c++){    
    cout << "Chain : "<<_sys.getChain(c).getChainId()<<endl;
    CartesianPoint thisChainGC(0.0,0.0,0.0);

    for (uint cterm=1; cterm <= 6; cterm++){
      Position &pos1 = _sys.getChain(c).getPosition(_sys.getChain(c).positionSize()-cterm);
      if (!pos1.atomExists("CA")) continue;
      cout << "ADDING POSITION: "<<pos1.getPositionId()<<endl;
      thisChainGC += pos1.getAtom("CA").getCoor();
    }
    cout << "Here"<<endl;
    thisChainGC /= 6;

    points.push_back(thisChainGC);
  }
    cout << "Here2"<<endl;
    // Compute a frame please!
  Frame aFrame;
  aFrame.computeFrameFrom3Points(points[0],points[1],points[2],true);

  return (aFrame);
}
