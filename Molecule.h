#ifndef MOLECULE_H
#define MOLECULE_H

#include <string>
#include <map>
#include <cctype>
using namespace std;

class Molecule{
    public:
        Molecule();
        void setName(string n);
        void setMolecularWeight(double mw);
        void setMass(double m);
        void setOccurence(int o);

        string getName();
        double getMolecularWeight();
        double getMass();
        int getOccurence();

    private:
        string name;
        double molecularWeight;
        double mass;
        int occurence;
};

#endif