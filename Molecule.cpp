#include "Molecule.h"

Molecule::Molecule(){
    name = "N";
    molecularWeight = 0;
    mass = 0;
    occurence = 1;
}


void Molecule::setName(string n){
    name = n;
}

void Molecule::setMolecularWeight(double mw){
    molecularWeight = mw;
}

void Molecule::setMass(double m){
    mass = m;
}

void Molecule::setOccurence(int o){
    occurence = o;
}

string Molecule::getName(){
    return name;
}

double Molecule::getMolecularWeight(){
    return molecularWeight;
}

double Molecule::getMass(){
    return mass;
}

int Molecule::getOccurence(){
    return occurence;
}