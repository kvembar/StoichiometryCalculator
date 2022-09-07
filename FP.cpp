/*
 *
 * TODOS/Caveats to fix:
 * #3) Each element in each compound appears once and only once.
 * #6) No checking for nonexistent elements.
 * Note: Doesn't check for balanced eq.
 * 
*/

#include "Molecule.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

//Calculates Molar Mass of a given formula. N is used for recursion purposes.

double MM(string formula, int N){
    //Adds one to the end of the formula so that the fxn may terminate properly.
    if(!isdigit(formula.back())){
        formula += '1';
    }

    //Reads in the atomic weights in the atomweight.txt file.
    ifstream fin("atomweight.txt");
    if (fin.fail()){
        cout << "File failed to open. Terminating..." << endl;
        return -1.0;
    }

    //Input of all element molecular masses and assignment within the hashmap.
    map<string, double> weight;
    string atom;
    double molar_mass;
    while(fin.good()){
        fin >> atom >> molar_mass;
        weight[atom] = molar_mass;
    }
    fin.close();

    //Calculation of molar mass of the element.
    double molarmass = 0; //What will be returned.
    int ind = 0; //holds current index of the string.
    string temp = ""; //temporary string to hold the element.

    while(ind < formula.length()){
        //TODO: Issue #3, Issue #6.

        //If you encounter an uppercase letter, its the start of a new element. 
        //At the very beginning of a compound or just after a number, 
        //it just adds a 0 and then assigns temp to the character there, so the edge case is handled.
        //If its in the middle, then there is only one of the element in question, so its added with any recursive multipliers handled.
        if(isupper(formula[ind])){
            molarmass += weight[temp]*N;
            temp = formula[ind];
        }
        //If you encounter a lowercase letter, its just another part of the element. Add it to temp and move on.
        else if(islower(formula[ind])){
            temp += formula[ind];
        }
        //If there is a digit, end of element and occurence number. Add it to molarmass.
        else if(isdigit(formula[ind])){
            molarmass += weight[temp] * (int)(formula[ind] - '0') * N;
            temp = ""; //Reset temp to the empty string.
        }
        //Recursive case with '(' and ')'.
        else if(formula[ind] == '('){

            if(temp != ""){
                molarmass += weight[temp] * N;
                temp = "";
                //If the string isn't empty, then we have some element we haven't dealt with yet.
                //Identical to uppercase case.
            }

            int layers = 1; //Keeps track of what layer parentheses we're on.

            ind++; //Moves us past the "("
            //While we are in the parentheses, add everything into temp and watch for the end of the parentheses.
            while((layers != 0) && (ind < formula.length())){
                temp += formula[ind];
                ind++;
                if(formula[ind] == '('){
                    layers++;
                }
                else if(formula[ind] == ')'){
                    layers--;
                }
            }
            ind++; //Go past the end ')'

            //Error handling of parentheses.
            if(layers != 0){
                return -1;
            }
            
             //Calculate MM of everything in parentheses and add it, handling case with and without number.
            if(isdigit(formula[ind])){
                molarmass += MM(temp,(formula[ind] - '0'));
            }
            else{
                molarmass += MM(temp,1);
            }
            temp = "";
            ind++;
        }
        else{
            return -1; //Error handling of anything that isn't alphanum or opening parentheses.
        }

        ind++;
    }

    return molarmass;
    
}

int main(){
    //txtfile input and checking for errors.
    cout << "Welcome to the Stoichiometry Calculator!" << endl;
    string txtfile;

    ifstream fin("formula.txt");

    if(fin.fail()){
        cerr << "File could not be opened. Terminating..." << endl;
        return -1;
    }

    //Reactants and products input, with reactants and products being the formulae names and the coeff variants being their coefficients
    vector<string> reactants;
    vector<int> coeffReactants;
    vector<string> products;
    vector<int> coeffProducts;
    bool reactantMode = true;

    string in; //To store the input of the variable

    while(fin.good()){
        fin >> in;

        //Skip of the = and +s within the txt file. If an = is encountered, switch to adding products instead of reactants.
        if(in == "="){
            reactantMode = false;
            continue;
        }

        if(in == "+"){
            continue;
        }

        //Determining the size and number of each of the elements and setting of various aspects of the molecule.
        if(reactantMode){
            if(isdigit(in[0])){
                string num = "";

                for(char i:in){
                    if(isdigit(i)){
                        num += i; //Stores numbers until there is no number left.
                        continue;
                    }
                    break;
                }

                coeffReactants.push_back(stoi(num));
                reactants.push_back( in.substr(num.length(), in.length() - num.length() ) );
            }
            else{
                coeffReactants.push_back(1);
                reactants.push_back( in );
            }
        }
        //Begin adding products now.
        else{
            if(isdigit(in[0])){
                string num = "";

                for(char i:in){
                    if(isdigit(i)){
                        num += i;
                        continue;
                    }
                    break;
                }

                coeffProducts.push_back(stoi(num));
                products.push_back( in.substr(num.length(), in.length() - num.length() ) );
            }
            else{
                coeffProducts.push_back(1);
                products.push_back( in );
            }
        }
        
    }
    fin.close();

    //Error checking for inclusion of =.
    if(reactants.size() == 0 || products.size() == 0){
        cerr << "Error occurred: Products or Reactants empty. Please check formatting." << endl;
        return -1;
    }

    //Formatting to "Molecule" format.
    vector<Molecule> molReactant(reactants.size());
    vector<Molecule> molProduct(products.size());
    double temp;

    //Setting of the various variables in the actual molecule format.
    cout << "Note: If the amount of material is in excess, enter in -1" << endl;
    for(int i = 0; i < reactants.size(); i++){
        molReactant.at(i).setName(reactants.at(i));
        molReactant.at(i).setOccurence(coeffReactants.at(i));

        temp = MM(reactants.at(i), 1);

        //Error checking in action.
        if(temp == -1){
            cerr << "Invalid chemical formula (invalid char detected): " << reactants.at(i) << endl;
            return -1;
        }

        molReactant.at(i).setMolecularWeight(temp);
        cout << "Please enter in the mass of the compound " << reactants.at(i) << "(g): ";

        cin >> temp;
        molReactant.at(i).setMass(temp);
    }

    for(int i = 0; i < products.size(); i++){  
        molProduct.at(i).setName(products.at(i));
        molProduct.at(i).setOccurence(coeffProducts.at(i));

        temp = MM(products.at(i), 1);

        if(temp == -1){
            cerr << "Invalid chemical formula (invalid char detected): " << reactants.at(i) << endl;
            return -1;
        }

        molProduct.at(i).setMolecularWeight(temp);
    }

    //Determining the Limiting Reactant based on how much product it produces. Stoich.
    string LR = "";
    int where;
    double productions = -1;

    for(int i = 0; i < reactants.size(); i++){
        //Skips any excess reactant specified by the user.
        if(molReactant.at(i).getMass() == -1){
            continue;
        }

        Molecule R = molReactant.at(i);
        Molecule P = molProduct.at(0);

        temp = (R.getMass())/(R.getMolecularWeight()) * (P.getOccurence())/(R.getOccurence());

        //Whichever produces the least is LR. productions == -1 handles first reactant by setting it as temp LR
        if((productions == -1) || (temp < productions)){
            LR = R.getName(); //Set LR to name
            where = i; //Keep where it is for output.
            productions = temp; //New min is this variable.
        }
    }

    //Use LR as basis to calculate the produced amount of the compounds in grams.
    Molecule basis = molReactant.at(where);
    double actualMass;
    for(int i = 0; i < reactants.size(); i++){
        if(i == where){
            continue;
        }
        //Same as before, but now we calculate EVERY SINGLE compound, reactant or product, using LR as the basis.
        Molecule P = molReactant.at(i);
        actualMass = (basis.getMass())/(basis.getMolecularWeight()) * (P.getOccurence())/(basis.getOccurence()) * (P.getMolecularWeight());
        molReactant.at(i).setMass(actualMass);
    }
    for(int i = 0; i < products.size(); i++){
        Molecule P = molProduct.at(i);
        actualMass = (basis.getMass())/(basis.getMolecularWeight()) * (P.getOccurence())/(basis.getOccurence()) * (P.getMolecularWeight());
        molProduct.at(i).setMass(actualMass);
    }

    //Output to text file called "results.txt."

    ofstream fout("results.txt");
    if(fout.fail()){
        cerr << "File could not be opened. Terminating..." << endl;
        return -1;
    }

    //Outputting to file and labeling limiting reactant
    //Also formatting it according to the longest element entered.
    fout << "Results: " << endl;
    fout << endl << "Reactants: " << endl;
    fout << setprecision(3) << fixed;

    //Max string calculated.
    int max_space = reactants.at(0).length();
    for(string i: reactants){
        if(i.length() > max_space) max_space = i.length();
    }
    for(string i: products){
        if(i.length() > max_space) max_space = i.length();
    }

    //Outputting to the file.
    for(int i = 0; i < reactants.size(); i++){

        fout << "Compound: " << left << setw(max_space) << reactants.at(i) 
        << " | Used Amount: " << molReactant.at(i).getMass() << " g / " 
        << molReactant.at(i).getMass() / MM(molReactant.at(i).getName(),1) << " mol";

        if(i == where){
            fout << " (LR)" << endl; //Label LR based on the 'where' variable in the LR calculating section.
        }
        else{
            fout << endl;
        }
    }

    fout << endl << "Products: " << endl;
    for(int i = 0; i < products.size(); i++){
        fout << "Compound: " << left << setw(max_space) << products.at(i) 
        << " | Produced Amount: " << molProduct.at(i).getMass() << " g / "
        << molProduct.at(i).getMass() / MM(molProduct.at(i).getName(),1) << " mol" << endl;
    }
    
    return 0; //Terminate program.
}