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
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  G4String name, symbol;             // a=mass of a mole;
  G4double a, z, density,abundance;            // z=mean number of protons;
  G4int iz, n;                       // iz=nb of protons  in an isotope;
                                    // n=nb of nucleons in an isotope;
  G4int ncomponents, natoms;

  G4UnitDefinition::BuildUnitsTable();

  // define Elements
  a = 1.01*g/mole;  G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);
  a = 12*g/mole;  G4Element* elC  = new G4Element(name="Carbon",symbol="C" , z= 6., a);
  a = 18*g/mole;  G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);
  a = 26.98*g/mole;  G4Element* elAl  = new G4Element(name="Aluminium"  ,symbol="Al" , z= 13., a);
  a = 208.98*g/mole;  G4Element* elBi  = new G4Element(name="Bismuth"  ,symbol="Bi" , z= 83., a);
  a = 72.64*g/mole;  G4Element* elGe  = new G4Element(name="Germanium"  ,symbol="Ge" , z= 32., a);
  // define an Element from isotopes, by relative abundance
  G4Isotope* Co59 = new G4Isotope(name="Co59", iz=27, n=59, a=58.9*g/mole);
  G4Isotope* Co60 = new G4Isotope(name="Co60", iz=27, n=60, a=60*g/mole);
  G4Element* elCo  = new G4Element(name="Cobalt", symbol="Co", ncomponents=2);
  elCo->AddIsotope(Co59, abundance= 90.*perCent);
  elCo->AddIsotope(Co60, abundance= 10.*perCent);

  // define simple materials
  density = 2.700*g/cm3;
  a = 26.98*g/mole;
  G4Material* Al = new G4Material(name="Aluminum", z=13., a, density);

  density = 7.13*g/cm3;
  G4Material* BGO= new G4Material(name="BGO", density, ncomponents=3);
  BGO->AddElement(elBi , natoms=4);
  BGO->AddElement(elGe, natoms=3);
  BGO->AddElement(elO , natoms=12);

  density = 6.44*g/cm3; // CoO
  G4Material* CoO= new G4Material(name="CoO", density, ncomponents=2);
  CoO->AddElement(elCo , natoms=1);
  CoO->AddElement(elO , natoms=1);

  density = 0.9*g/cm3; // 塑料
  G4Material* Plastic= new G4Material(name="Plastic", density, ncomponents=2);
  Plastic->AddElement(elC , natoms=1);
  Plastic->AddElement(elH , natoms=1);

  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;

  //
  // World
  //
  G4double world_sizeXY = 5*1.618*cm;
  G4double world_sizeZ  = 5*cm;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  G4Box* solidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

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
  // 重新定义Envelope用作放射源基材，即射线源从中抽样

  G4double env_sizeXY = 1.5*cm, env_sizeZ = 0.2*cm;
  G4Material* env_mat = Plastic; 
  // G4Material* env_mat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL"); 
  env_pos = G4ThreeVector(0,0,0.5*env_sizeZ);

  // G4Material* env_mat = CoO;
  G4Box* solidEnv =
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size

  G4LogicalVolume* logicEnv =
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name

  new G4PVPlacement(0,                       //no rotation
                    env_pos,         //at (0,0,0.5*env_sizeZ)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //
  // Al Frame
  //
  G4Material* Alframe_mat = Al;
  G4double Alframe_sizeX = 2.7*cm, Alframe_sizeY = 6*cm, Alframe_sizeZ = 1*cm;
  G4double AlframeThickness = 0.1*cm;
  G4ThreeVector Alframe_pos = G4ThreeVector(0, 0, -0.5*Alframe_sizeZ);
  G4Box* solidAlframe = new G4Box("Al_Frame", 0.5*Alframe_sizeX, 0.5*Alframe_sizeY, 0.5*Alframe_sizeZ); 
  G4LogicalVolume* logicAlframe = new G4LogicalVolume(solidAlframe, Alframe_mat,"Alframe");
  new G4PVPlacement(0,Alframe_pos,logicAlframe,"Alframe",logicWorld,false,0,checkOverlaps);

  //
  // BGO
  //
  G4Material* BGO_mat = BGO;
  G4double BGO_sizeX = Alframe_sizeX-2*AlframeThickness;
  G4double BGO_sizeY = Alframe_sizeY-2*AlframeThickness;
  G4double BGO_sizeZ = Alframe_sizeZ-2*AlframeThickness;
  G4ThreeVector BGO_pos = G4ThreeVector(0, 0, 0); // 相对母体的位置
  G4Box* solidBGO = new G4Box("BGO_Scintillator", 0.5*BGO_sizeX, 0.5*BGO_sizeY, 0.5*BGO_sizeZ); 
  G4LogicalVolume* logicBGO = new G4LogicalVolume(solidBGO, BGO_mat,"BGO");
  new G4PVPlacement(0,BGO_pos,logicBGO,"BGO",logicAlframe,false,0,checkOverlaps);

  // Set BGO as scoring volume
  //
  fScoringVolume = logicBGO;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
