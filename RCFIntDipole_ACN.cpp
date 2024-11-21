#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

struct Atom {
    std::string symbol;
    double x, y, z;
};

struct Vector {
    double x, y, z;
};

struct Molecule {
    Atom M;
    Atom C;
    Atom N;
    Vector CM;
    Vector CN;
    Vector Mu;
};

std::vector<std::vector<Atom>> read_frames(const std::string& filename, int num_frames) {
    std::ifstream file(filename);
    int num_atoms = 1152;
    double dum;
    std::string dummy_line;
    std::vector<std::vector<Atom>> all_frames;

    for (int i = 0; i < num_frames; ++i) {
//        file >> num_atoms;
//        std::getline(file, dummy_line);  // read the rest of the first line
//        std::getline(file, dummy_line);  // read the second line into a dummy variable

        std::vector<Atom> atoms(num_atoms);
        for (Atom& atom : atoms) {
            file >> atom.symbol >> atom.x >> atom.y >> atom.z >> dum >> dum >> dum;
        }

        // Now atoms contains the atom symbols and coordinates for this frame
        //Add this frame to all_frames
        all_frames.push_back(atoms);

	if (i%10000 == 0) {
		std::cout << i << "Frames Read" << std::endl;
	}

        }
        //Now all_frames contains all the frames from the file
        return all_frames;
    }
        
void write_frames(const std::string& filename, const std::vector<std::vector<Atom>>& all_frames) {
    int frame_index = 0;    
    std::ofstream file(filename);
    for (const auto& frame : all_frames) {
        file << frame.size() << "\n";
        file << " Atoms. Timestep: " << frame_index*10000 << "\n";  // You can replace this with your own comment
        for (const auto& atom : frame) {
            file << atom.symbol << " " << std::setprecision(8) << atom.x << " " << atom.y << " " << atom.z << "\n";
        }
        ++frame_index;
    }
}

Vector compute_normalized_vector(const Atom& atom1, const Atom& atom2) {
    Vector vec = {atom2.x - atom1.x, atom2.y - atom1.y, atom2.z - atom1.z};
    double magnitude = std::sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
    return {vec.x / magnitude, vec.y / magnitude, vec.z / magnitude};
}

Vector compute_simple_vector(const Atom& atom1, const Atom& atom2) {
    Vector vec = {atom2.x - atom1.x, atom2.y - atom1.y, atom2.z - atom1.z};
    return vec;
}

double compute_norm_vector(Vector vec) {
    double magnitude = std::sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
    return magnitude;
}


Vector compute_dipole_vector(const Molecule& molecule) {
    // Input charges
    double q_M = 0.4238;
    double q_C = 0.4238;
    double q_N = 0.4238;

    // Compute the dipole vector in e·Å using CM and CN vectors
    Vector dipole = {
        q_M * molecule.CM.x +  q_N * molecule.CN.x,
        q_M * molecule.CM.y +  q_N * molecule.CN.y,
        q_M * molecule.CM.z +  q_N * molecule.CN.z
    };

    double magnitude = std::sqrt(dipole.x * dipole.x + dipole.y * dipole.y + dipole.z * dipole.z); // Do I need to do this?
    return {dipole.x / magnitude, dipole.y / magnitude, dipole.z / magnitude};
}

std::vector<std::vector<Molecule>> compute_molecules(const std::vector<std::vector<Atom>>& all_frames) {
    std::vector<std::vector<Molecule>> molecules_data(all_frames.size());  // Intialization of a empty vector with type molecule 

    // Iterate over frames
    for (size_t t = 0; t < all_frames.size(); ++t) {
        // Iterate over atoms (assuming each acn molecule is represented by 3 atoms: M, C, N)
        for (size_t i = 0; i < all_frames[t].size(); i += 3) {
            Molecule molecule = {all_frames[t][i], all_frames[t][i + 1], all_frames[t][i + 2]};
            // compute CM and CM
            molecule.CM = compute_simple_vector(molecule.C, molecule.M);
            molecule.CN = compute_simple_vector(molecule.C, molecule.N);
            // compute dipole
	    molecule.Mu = compute_dipole_vector(molecule);
	    // std::cout << compute_norm_vector(molecule.Mu) << std::endl;
            // Append the molecule to molecules_in_slab
            molecules_data[t].push_back(molecule);
        }
    }

    return molecules_data;
}

double dot_product(const Vector& v1, const Vector& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

std::vector<long double> compute_tcf_dipole(const std::vector<std::vector<Molecule>>& molecules_data, int nframes, int ngap) {
    std::vector<long double> corr(nframes, 0.0);
    std::vector<int> num(nframes, 0);
    for (size_t j = 0; j < molecules_data[0].size(); j++) {  // j is over molecules
        for (int i = 0; i < nframes; i += ngap) {    // T0    // i is over frames
            int kmax = nframes;    // Tmax
            int k = i;
            while ((k < kmax && molecules_data[i][j].C.z >= 0.0 && molecules_data[k][j].C.z >= 0.0 && molecules_data[i][j].C.z <= 5.14 && molecules_data[k][j].C.z <= 5.14) || (k < kmax && molecules_data[i][j].C.z >= 37.90 && molecules_data[k][j].C.z >= 37.90 && molecules_data[i][j].C.z <= 45.0 && molecules_data[k][j].C.z <= 45.0) )             {
                corr[k-i] += dot_product(molecules_data[i][j].Mu, molecules_data[k][j].Mu);
                num[k-i]+=1;
                //std::cout << "yes" << std::endl ;
                k++;
            }
         //break;
        }
      // break;
    }

    // Average over all molecules
    for (int i = 0; i < nframes; i++) {
        if (num[i] != 0) {
            corr[i] /= num[i];
        }
    }
   
    // Normalization by T0 
    //for (size_t i = 0; i < corr.size(); ++i) {
    //    corr[i] = corr[i]/corr[0];
    //}

    return corr;
}


void write_to_file(const std::vector<long double>& corr1, const std::string& filename) {
    std::ofstream file(filename);

    if (file.is_open()) {
        for (size_t i = 0; i < corr1.size(); ++i) {
            file << i << " " << std::setprecision(14) << corr1[i] << "\n";
        }
        file.close();
    } else {
        std::cout << "Unable to open file: " << filename << std::endl;
    }
}



int main() {
    std::string filename = "../LongFreqSavedRun/posvel.xyz";
    int num_frames = 300000;  // replace with the actual number of frames you want to read
        
    std::vector<std::vector<Atom>> all_frames = read_frames(filename, num_frames);

    // std::string output_filename = "output.xyz";  // replace with your actual output filename
    // write_frames(output_filename, all_frames);    
   
    std::vector<std::vector<Molecule>> molecules_MCN = compute_molecules(all_frames); 
    std::vector<long double> tcf = compute_tcf_dipole(molecules_MCN, num_frames, 2);
    write_to_file(tcf, "Intdipole.dat");
    return 0;
    }
