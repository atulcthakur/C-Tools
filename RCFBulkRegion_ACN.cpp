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
    Vector MN;
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

std::vector<std::vector<Molecule>> compute_molecules(const std::vector<std::vector<Atom>>& all_frames) {
    std::vector<std::vector<Molecule>> molecules_data(all_frames.size());

    // Iterate over frames
    for (size_t t = 0; t < all_frames.size(); ++t) {
        // Iterate over atoms (assuming each molecule is represented by 3 atoms: M, C, N)
        for (size_t i = 0; i < all_frames[t].size(); i += 3) {
            Molecule molecule = {all_frames[t][i], all_frames[t][i + 1], all_frames[t][i + 2]};
            // compute OH1 and OH2
            molecule.MN = compute_normalized_vector(molecule.M, molecule.N);
            // Append the molecule to molecules_in_slab
            molecules_data[t].push_back(molecule);
        }
    }

    return molecules_data;
}

double dot_product(const Vector& v1, const Vector& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}


std::vector<long double> compute_time_correlation(const std::vector<std::vector<Molecule>>& molecules_data, int nframes, int ngap) {
    std::vector<long double> corr(nframes, 0.0);
    std::vector<int> num(nframes, 0);
    for (size_t j = 0; j < molecules_data[0].size(); j++) {  // j is over molecules
	std::cout << "Doing Carbon " << j << std::endl ;
        for (int i = 0; i < nframes; i += ngap) {    // T0    // i is over frames
            int kmax = nframes;    // Tmax
            int k = i;
            while ((k < kmax && molecules_data[i][j].C.z >= 18.53 && molecules_data[k][j].C.z >= 18.53 && molecules_data[i][j].C.z <= 21.53 && molecules_data[k][j].C.z <= 21.53) || (k < kmax && molecules_data[i][j].C.z >= 21.53 && molecules_data[k][j].C.z >= 21.53 && molecules_data[i][j].C.z <= 24.53 && molecules_data[k][j].C.z <= 24.53) )             {
                corr[k-i] += dot_product(molecules_data[i][j].MN, molecules_data[k][j].MN);
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

    return corr;
}

std::vector<long double> compute_P2time_correlation(const std::vector<std::vector<Molecule>>& molecules_data, int nframes, int ngap) {
    std::vector<long double> corr(nframes, 0.0);
    std::vector<int> num(nframes, 0);
    for (size_t j = 0; j < molecules_data[0].size(); j++) {  // j is over molecules
        std::cout << "Doing Carbon " << j << std::endl ;
        for (int i = 0; i < nframes; i += ngap) {    // T0    // i is over frames
            int kmax = nframes;    // Tmax
            int k = i;
            while ((k < kmax && molecules_data[i][j].C.z >= 18.53 && molecules_data[k][j].C.z >= 18.53 && molecules_data[i][j].C.z <= 21.53 && molecules_data[k][j].C.z <= 21.53) || (k < kmax && molecules_data[i][j].C.z >= 21.53 && molecules_data[k][j].C.z >= 21.53 && molecules_data[i][j].C.z <= 24.53 && molecules_data[k][j].C.z <= 24.53) )             {
                double MN_dot = dot_product(molecules_data[i][j].MN, molecules_data[k][j].MN);
                corr[k-i] += (1.5 * MN_dot*MN_dot - 0.5);
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
    // int num_frames = 1000000;  // replace with the actual number of frames you want to read
    int num_frames = 300000;  // replace with the actual number of frames you want to read
        
    std::vector<std::vector<Atom>> all_frames = read_frames(filename, num_frames);

    // std::string output_filename = "output.xyz";  // replace with your actual output filename
    // write_frames(output_filename, all_frames);    
   
    std::vector<std::vector<Molecule>> molecules_MN = compute_molecules(all_frames); 
    std::vector<long double> tcf = compute_time_correlation(molecules_MN, num_frames, 2);
    write_to_file(tcf, "BulkRegion.dat");
    
    std::vector<long double> P2tcf = compute_P2time_correlation(molecules_MN, num_frames, 2);
    write_to_file(P2tcf, "P2BulkRegion.dat");

    return 0;
    }
