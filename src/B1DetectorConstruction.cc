//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B1DetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Tubs.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  G4double env_sizeXY = 20*cm, env_sizeZ = 30*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;


  //G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  G4double a, z, density, temperature, pressure;
  G4String name, symbol;    
  

  G4int ncomponents, natoms;
  // Define Boron 
  G4Element* elementB=new G4Element(name="Boron", symbol="B", z=5., a=10.811*g/mole);

  // Define Carbon
  G4Element* elementC=new G4Element(name="Carbon", symbol="C", z=6., a=12.011*g/mole);

  // Define Oxygen 
  G4Element* elementO=new G4Element(name="Oxygen", symbol="O", z=8., a=15.9994*g/mole);

  // Define Sodium
  G4Element* elementNa = new G4Element(name="Sodium", symbol="Na", z=11., a=22.9898*g/mole);

  // Define Silicon
  G4Element* elementSi = new G4Element(name="Silicon", symbol="Si", z=14., a=28.09*g/mole);

  // Define Iron
  G4Element* elementFe = new G4Element(name="Iron", symbol="Fe", z=26., a=55.845*g/mole);

  // Define Copper
  G4Element* elementCu = new G4Element(name="Copper", symbol="Cu", z=29., a=63.5460*g/mole);

  // Define Germanium   
  G4Element* elementGe = new G4Element(name="Germanium", symbol="Ge", z=32., a=72.61*g/mole);

  // Define Iodine   
  G4Element* elementI = new G4Element(name="Iodine", symbol="I", z=53., a=126.90447*g/mole);

  // Define Lead
  G4Element* elementPb = new G4Element(name="Lead",symbol="Pb", z=82., a=207.2*g/mole);
  
  // Define Phosphorus
  G4Element* elementP = new G4Element(name="Phosphorus",symbol="P", z=15., a=30.97*g/mole);
  
  // Define Sulfur
  G4Element* elementS = new G4Element(name="Sulfur",symbol="S", z=16., a=32.065*g/mole);
  
  // Define Nickel
  G4Element* elementNi = new G4Element(name="Nickel",symbol="Ni", z=28., a=58.6934*g/mole);

  // Define Chromium
  G4Element* elementCr = new G4Element(name="Chromium",symbol="Cr", z=24., a=51.9961*g/mole);

  // Define Manganese
  G4Element* elementMn = new G4Element(name="Manganese",symbol="Mn", z=25., a=54.9380*g/mole);

  // Define Zinc
  G4Element* elementZn = new G4Element(name="Zinc",symbol="Zn", z=30., a=65.38*g/mole);

  //     
  // World
  //
  //G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeX = 347.94*cm;
  G4double world_sizeY = 403.94*cm;
  G4double world_sizeZ  = 1100.0*cm;
  //G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Envelope
  //  
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
  G4VPhysicalVolume* physicalWorld = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 
  logicEnv->SetVisAttributes(G4VisAttributes::GetInvisible());

  G4double Tower_nZcut = 2;
  G4double Tower_zPcut[] = {-126.231,  126.231};
  G4double Tower_rIcut[] = {0, 0};
  G4double Tower_rOcut[] = {57.34, 42.45};
  

  // Vacuum
  density     = universe_mean_density;   //from PhysicalConstants.h
  pressure    = 1.0e-19*pascal; 
  temperature = 2.73*kelvin;
  G4Material* Vacuum = new G4Material("Vacuum", z=1.0, a=1.01*g/mole, density,
              kStateGas, temperature, pressure);

  G4ThreeVector pos1 = G4ThreeVector(0, 0, -6.5*cm);
  
  G4RotationMatrix r180Rotation;    // flip towers over
  r180Rotation.rotateY(180.*deg);
  G4Transform3D towerflip(r180Rotation, pos1);

  G4Polyhedra* solidTower1 = new G4Polyhedra("Tower1_S",0.*deg,360.*deg,6,
    Tower_nZcut,Tower_zPcut,Tower_rIcut,Tower_rOcut);

  G4LogicalVolume* logicShape1 =                         
    new G4LogicalVolume(solidTower1,         //its solid
                        Vacuum,          //its material
                        "Shape1");           //its name
               
  // new G4PVPlacement(0,                       //no rotation
  //                   pos1,                    //at position
  //                   logicShape1,             //its logical volume
  //                   "Shape1",                //its name
  //                   logicEnv,                //its mother  volume
  //                   false,                   //no boolean operation
  //                   0,                       //copy number
  //                   checkOverlaps);          //overlaps checking

  G4VPhysicalVolume* physicalTower1 = new G4PVPlacement(towerflip, // flipping with rotation here kills it, even putting no rotation here kills it.  Must use a G4Transform3D, but it works perfectly.
                /*positionTower,*/           // in towerflip
                "Shape1_P",             // name
                logicShape1,         // the logical volume to use
                physicalWorld,  // the mother volume
                false,               // no boolean operation
                0); 

G4VisAttributes* VisAtttmp = new G4VisAttributes(G4Colour(0.4,0.4,0.4,0.6));
//logicShape1->SetVisAttributes(G4VisAttributes::GetInvisible());
logicShape1->SetVisAttributes(VisAtttmp);



G4double zPcut[2];
G4double rIcut[2];
G4double rOcut[2];

G4double zPclc[2];
G4double rIclc[2];
G4double rOclc[2];


G4double rIctu;
G4double rOctu;
G4double zHctu;

G4double rIrb;
G4double rOrb;
G4double zHrb;

G4double zPudh[2];
G4double rIudh[2];
G4double rOudh[2];

G4double zPsdh[2];
G4double rIsdh[2];
G4double rOsdh[2];

G4double zPldh[2];
G4double rIldh[2];
G4double rOldh[2];

G4double Zip_z;

// "Copper Upper Tower" Parameters
zPcut[0]=-5.3108*cm;  zPcut[1]= 5.3108*cm;
rIcut[0]=3.0771*cm;   rIcut[1]=3.0771*cm;
rOcut[0]=3.9428*cm;   rOcut[1]=3.9428*cm;

// Lower cap of the "copper upper tower" Parameters
zPclc[0]=-0.2425*cm; zPclc[1]= 0.2425*cm;
rIclc[0]=0*cm;       rIclc[1]=0*cm;
rOclc[0]=3.9428*cm;  rOclc[1]=3.9428*cm;

// Connector Tube Parameters
rIctu = 2.6924*cm;   // Rinner of Cu tube
rOctu = 2.794*cm;   // Router of Cu tube
zHctu = 4.48691*cm; // Half Length of Cu tube

// Base of Connector Tube Parameters
rIrb = 2.8067*cm;   // Rinner of Cu tube
rOrb = 3.429*cm;   // Router of Cu tube
zHrb = 0.2405*cm; // Half Length of Cu tube

// Upper cap of the detector housing Parameters
zPudh[0]=-0.3277*cm; zPudh[1]= 0.3277*cm;
rIudh[0]=0*cm;       rIudh[1]=0*cm;
rOudh[0]=5.334*cm; rOudh[1]=5.334*cm;

// The side detector housing
zPsdh[0]=-1.8174*cm;  zPsdh[1]= 1.8174*cm;
rIsdh[0]= 5.1308*cm;  rIsdh[1]= 5.1308*cm;
rOsdh[0]= 5.334*cm; rOsdh[1]= 5.334*cm;

// The lower cap of the detector housing
zPldh[0]=-0.3175*cm;  zPldh[1]= 0.3175*cm;
rIldh[0]=0*cm;         rIldh[1]=0*cm;
rOldh[0]=5.334*cm;   rOldh[1]=5.334*cm;

Zip_z = 1.312*2.54*cm; //Corrected from 2.54*cm; 

// Visualization attributes
  G4VisAttributes* VisAttCu1 = new G4VisAttributes(G4Colour(180/255.,140/255.,60/255.)); VisAttCu1->SetForceSolid(true); 
  G4VisAttributes* VisAttCu2 = new G4VisAttributes(G4Colour(220/255.,165/255.,55/255.)); VisAttCu2->SetForceSolid(true); 
  G4VisAttributes* VisAttCu3 = new G4VisAttributes(G4Colour(210/255.,165/255.,70/255.)); VisAttCu3->SetForceSolid(true); 
  G4VisAttributes* VisAttCu4 = new G4VisAttributes(G4Colour(128/255.,64/255.,0/255.)); VisAttCu4->SetForceSolid(true); 


// Copper
  G4Material* Copper = new G4Material(name="Copper", density = 8.920*g/cm3, ncomponents=1);
  Copper->AddElement(elementCu, natoms=1);
  
  G4double nZcut = 2;
  G4double nZclc = 2;
  G4double nZudh = 2;
  G4double nZsdh = 2;
  G4double nZldh = 2;

  //-------------------------------------------------------------------------
  //Working from the top down: this is the `copper upper tower', mass 2.0 kg
  //-------------------------------------------------------------------------
  G4ThreeVector position_top = G4ThreeVector(0,0,+Tower_zPcut[1] - zPcut[1]);
  G4cout << "z shift for Tower Guts: " << Tower_zPcut[1] - zPcut[1] << G4endl;
  G4cout << "Height of copper upper tower: " << zPcut[1] - zPcut[0] << G4endl;
  //spandey
  G4cout << "nZcut of copper upper tower: " << nZcut << G4endl;

  G4Polyhedra* cu_tops=new G4Polyhedra("cu_top",0.*deg,360.*deg,6,nZcut,zPcut,rIcut,rOcut);
  G4LogicalVolume* cu_topl1 = new G4LogicalVolume(cu_tops,Copper,"cutl1");
  G4PVPlacement* cu_topp1 = new G4PVPlacement(0,position_top,"cutp1",cu_topl1,physicalTower1,false,0);
  cu_topl1->SetVisAttributes(VisAttCu1); 


  //----------------------------------------------------------------------
  //Next, the lower cap of the `copper upper tower', mass 0.233 kg/tower
  //----------------------------------------------------------------------
  G4ThreeVector position_clc = G4ThreeVector(0,0,+Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - zPclc[1]);
  G4cout << "z shift for Tower lower cap: " << Tower_zPcut[1] - (zPcut[1]-zPcut[0]) << G4endl;
  G4cout << "Height of copper lower cap: " << zPclc[1] - zPclc[0] << G4endl;
  //spandey
  G4cout << "nZclc of copper lower cap: " << nZclc << G4endl;
  G4Polyhedra* cu_clcs=new G4Polyhedra("cu_clc",0.*deg,360.*deg,6,nZclc,zPclc,rIclc,rOclc);
  G4LogicalVolume* cu_clcl1 = new G4LogicalVolume(cu_clcs,Copper,"clcl1");
  G4PVPlacement* cu_clcp1 = new G4PVPlacement(0,position_clc,"clcp1",cu_clcl1,physicalTower1,false,0);
  cu_clcl1->SetVisAttributes(VisAttCu2); 


  //------------------------------------------------------------------
  //Next, connector tube, mass 0.140 kg/tower (next ring excluded...)
  //------------------------------------------------------------------
  G4ThreeVector position_ctu = G4ThreeVector(0,0,+Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - zHctu);
  G4cout << "z shift for spool: " << Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - zHctu << G4endl;
  G4cout << "Height of spool: " << 2*zHctu << G4endl;
  G4Tubs* cu_ctus=new G4Tubs("cu_ctu",rIctu,rOctu,zHctu,0.*deg,360.*deg);
  G4LogicalVolume* cu_ctul1=new G4LogicalVolume(cu_ctus,Copper,"ctul1");
  G4PVPlacement* cu_ctup1=new G4PVPlacement(0,position_ctu,"ctup1",cu_ctul1,physicalTower1,false,0);
  cu_ctul1->SetVisAttributes(VisAttCu3); 

  //--------------------------------------------------------------
  //Next, ring at base of connector tube, mass 0.051 kg/tower
  //--------------------------------------------------------------
  G4ThreeVector position_rb = G4ThreeVector(0,0, Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu); // QUESTION : Is this ring centered on the bottom of the tower ??
  G4cout << "z shift for base ring: " << Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu << G4endl;
  G4cout << "Height of base ring: " << 2*zHrb << G4endl;
  G4Tubs* cu_rbs=new G4Tubs("cu_rb",rIrb,rOrb,zHrb,0.*deg,360.*deg);
  G4LogicalVolume* cu_rbl1=new G4LogicalVolume(cu_rbs,Copper,"crbl1");
  G4PVPlacement* cu_rbp1=new G4PVPlacement(0,position_rb,"crbp1",cu_rbl1,physicalTower1,false,0);
  cu_rbl1->SetVisAttributes(VisAttCu1); 

   //------------------------------------------------------------------
  //Next, the upper cap of the detector housing, mass 0.071 kg/tower
  // (about 0.044 kg/tower actually included in ring at base above)
  //------------------------------------------------------------------
  G4ThreeVector position_udh = G4ThreeVector(0,0, Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu - zHrb - zPudh[1]); 
  G4cout << "z shift for upper det housing cap: " << Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu - zHrb - zPudh[1] << G4endl;
  G4cout << "Height of upper det housing cap: " << zPudh[1] - zPudh[0] << G4endl;
  //spandey
  G4cout << "nZudh of upper det housing cap: " << nZudh << G4endl;
  G4Polyhedra* cu_udhs=new G4Polyhedra("cu_udh",0.*deg,360.*deg,6,nZudh,zPudh,rIudh,rOudh);
  G4LogicalVolume* cu_udhl1 = new G4LogicalVolume(cu_udhs,Copper,"udhl1");
  G4PVPlacement* cu_udhp1 = new G4PVPlacement(0, position_udh,"udhp1",cu_udhl1,physicalTower1,false,0);
  cu_udhl1->SetVisAttributes(VisAttCu2);


  //--------------------------------------------------------------
  //Next, the side detector housing, mass 0.384 kg/tower
  //--------------------------------------------------------------
  G4ThreeVector position_sdh = G4ThreeVector(0,0, Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu - zHrb - (zPudh[1]-zPudh[0]) - zPsdh[1]); 
  G4cout << "z shift for side det housing: " << Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu - zHrb - (zPudh[1]-zPudh[0]) - zPsdh[1] << G4endl;
  G4cout << "Height of side det housing: " << zPsdh[1] - zPsdh[0] << G4endl;
  //spandey
  G4cout << "nZsdh of side det housing: " << nZsdh << G4endl;
  G4Polyhedra* cu_sdhs=new G4Polyhedra("cu_sdh",0.*deg,360.*deg,6,nZsdh,zPsdh,rIsdh,rOsdh);
  G4LogicalVolume* cu_sdhl1 = new G4LogicalVolume(cu_sdhs,Copper,"sdhl1");
  G4PVPlacement* cu_sdhp1 = new G4PVPlacement(0,position_sdh,"sdhp1",cu_sdhl1,physicalTower1,false,0);
  cu_sdhl1->SetVisAttributes(VisAttCu3);

  //------------------------------------------------------------------
  //Next, the lower cap of the detector housing, mass 0.078 kg/tower
  //------------------------------------------------------------------
  G4ThreeVector position_ldh = G4ThreeVector(0,0, Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu - zHrb - (zPudh[1]-zPudh[0]) - (zPsdh[1]-zPsdh[0]) - zPldh[1] ); 
  G4cout << "z shift for lower det housing cap: " << Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu - zHrb - (zPudh[1]-zPudh[0]) - (zPsdh[1]-zPsdh[0]) - zPldh[1] << G4endl;
  G4cout << "Height of lower det housing cap: " << zPldh[1] - zPldh[0] << G4endl;
  //spandey
  G4cout << "nZldh of lower det housing: " << nZldh << G4endl;
  G4Polyhedra* cu_ldhs=new G4Polyhedra("cu_ldh",0.*deg,360.*deg,6,nZldh,zPldh,rIldh,rOldh);
  G4LogicalVolume* cu_ldhl1 = new G4LogicalVolume(cu_ldhs,Copper,"ldhl1");
  G4PVPlacement* cu_ldhp1 = new G4PVPlacement(0,position_ldh,"ldhp1",cu_ldhl1,physicalTower1,false,0);
  cu_ldhl1->SetVisAttributes(VisAttCu2);
  //cu_sdhl1->SetVisAttributes(G4VisAttributes::Invisible);  // Make Invisible







  // Side Coax
  G4double sdcx_widthH = 1.4*cm;
  G4double sdcx_thicknessH = 0.065*cm;
  G4double sdcx_lenH[6];
  sdcx_lenH[0]=5.294*cm; //4.86*cm;
  sdcx_lenH[1]=5.294*cm; //5.64*cm;
  sdcx_lenH[2]=5.294*cm; //6.31*cm;
  sdcx_lenH[3]=5.294*cm; //6.99*cm;
  sdcx_lenH[4]=5.294*cm; //7.66*cm;
  sdcx_lenH[5]=5.294*cm; //8.34*cm;
  G4double sdcx_radius = 4.18*cm;

  //------------------------------------------------------------------
  //Next, the Side Coaxes
  //------------------------------------------------------------------

  //rotation matrices for side coaxes                                                     
  G4ThreeVector newz(0.0,0.0,1.0);  //z-axis doesn't change for any coax                  
  
  G4ThreeVector newx1(-cos(60*deg),sin(60*deg),0.0);
  G4ThreeVector newy1(-sin(60*deg),-cos(60*deg),0.0);
  G4RotationMatrix* coaxrotation1=new G4RotationMatrix(newx1,newy1,newz);
  
  G4ThreeVector newx2(-cos(60*deg),-sin(60*deg),0.0);
  G4ThreeVector newy2(sin(60*deg),-cos(60*deg),0.0);
  G4RotationMatrix* coaxrotation2=new G4RotationMatrix(newx2,newy2,newz);

  G4Box* solidcoax1=new G4Box("scoax1", sdcx_widthH, sdcx_thicknessH, sdcx_lenH[1-1]);
  G4Box* solidcoax2=new G4Box("scoax2", sdcx_widthH, sdcx_thicknessH, sdcx_lenH[2-1]);
  G4Box* solidcoax3=new G4Box("scoax3", sdcx_widthH, sdcx_thicknessH, sdcx_lenH[3-1]);
  G4Box* solidcoax4=new G4Box("scoax4", sdcx_widthH, sdcx_thicknessH, sdcx_lenH[4-1]);
  G4Box* solidcoax5=new G4Box("scoax5", sdcx_widthH, sdcx_thicknessH, sdcx_lenH[5-1]);
  G4Box* solidcoax6=new G4Box("scoax6", sdcx_widthH, sdcx_thicknessH, sdcx_lenH[6-1]);
  G4LogicalVolume* logiccoax1=new G4LogicalVolume(solidcoax1,Copper,"LogCoax1");
  G4LogicalVolume* logiccoax2=new G4LogicalVolume(solidcoax2,Copper,"LogCoax2");
  G4LogicalVolume* logiccoax3=new G4LogicalVolume(solidcoax3,Copper,"LogCoax3");
  G4LogicalVolume* logiccoax4=new G4LogicalVolume(solidcoax4,Copper,"LogCoax4");
  G4LogicalVolume* logiccoax5=new G4LogicalVolume(solidcoax5,Copper,"LogCoax5");
  G4LogicalVolume* logiccoax6=new G4LogicalVolume(solidcoax6,Copper,"LogCoax6");
  G4LogicalVolume* logiccoax7=new G4LogicalVolume(solidcoax1,Copper,"LogCoax7");
  G4LogicalVolume* logiccoax8=new G4LogicalVolume(solidcoax2,Copper,"LogCoax8");
  G4LogicalVolume* logiccoax9=new G4LogicalVolume(solidcoax3,Copper,"LogCoax9");
  G4LogicalVolume* logiccoax10=new G4LogicalVolume(solidcoax4,Copper,"LogCoax10");
  G4LogicalVolume* logiccoax11=new G4LogicalVolume(solidcoax5,Copper,"LogCoax11");
  G4LogicalVolume* logiccoax12=new G4LogicalVolume(solidcoax6,Copper,"LogCoax12");
  
  // Visualization line
  logiccoax1->SetVisAttributes(VisAttCu4);
  logiccoax2->SetVisAttributes(VisAttCu4);
  logiccoax3->SetVisAttributes(VisAttCu4);
  logiccoax4->SetVisAttributes(VisAttCu4);
  logiccoax5->SetVisAttributes(VisAttCu4);
  logiccoax6->SetVisAttributes(VisAttCu4);
  //     
  // Shape 1
  //  
  // G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  // G4ThreeVector pos1 = G4ThreeVector(0, 2*cm, -7*cm);
        
  // // Conical section shape       
  // G4double shape1_rmina =  0.*cm, shape1_rmaxa = 2.*cm;
  // G4double shape1_rminb =  0.*cm, shape1_rmaxb = 4.*cm;
  // G4double shape1_hz = 3.*cm;
  // G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
  // G4Cons* solidShape1 =    
  //   new G4Cons("Shape1", 
  //   shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb, shape1_hz,
  //   shape1_phimin, shape1_phimax);
                      
  // G4LogicalVolume* logicShape1 =                         
  //   new G4LogicalVolume(solidShape1,         //its solid
  //                       shape1_mat,          //its material
  //                       "Shape1");           //its name
               
  // new G4PVPlacement(0,                       //no rotation
  //                   pos1,                    //at position
  //                   logicShape1,             //its logical volume
  //                   "Shape1",                //its name
  //                   logicEnv,                //its mother  volume
  //                   false,                   //no boolean operation
  //                   0,                       //copy number
  //                   checkOverlaps);          //overlaps checking

  // //     
  // // Shape 2
  // //
  // G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
  // G4ThreeVector pos2 = G4ThreeVector(0, -1*cm, 7*cm);

  // // Trapezoid shape       
  // G4double shape2_dxa = 12*cm, shape2_dxb = 12*cm;
  // G4double shape2_dya = 10*cm, shape2_dyb = 16*cm;
  // G4double shape2_dz  = 6*cm;      
  // G4Trd* solidShape2 =    
  //   new G4Trd("Shape2",                      //its name
  //             0.5*shape2_dxa, 0.5*shape2_dxb, 
  //             0.5*shape2_dya, 0.5*shape2_dyb, 0.5*shape2_dz); //its size
                
  // G4LogicalVolume* logicShape2 =                         
  //   new G4LogicalVolume(solidShape2,         //its solid
  //                       shape2_mat,          //its material
  //                       "Shape2");           //its name
               
  // new G4PVPlacement(0,                       //no rotation
  //                   pos2,                    //at position
  //                   logicShape2,             //its logical volume
  //                   "Shape2",                //its name
  //                   logicEnv,                //its mother  volume
  //                   false,                   //no boolean operation
  //                   0,                       //copy number
  //                   checkOverlaps);          //overlaps checking
                
  // // Set Shape2 as scoring volume
  // //
  // fScoringVolume = logicShape2;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
