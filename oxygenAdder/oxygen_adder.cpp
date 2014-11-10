#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iterator>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

bool check_oxygens = false;

//
// A quick and dirty utility which adds oxygens to backbones of proteins.
// The utility uses PDB format.
//
// Usage forms: cat input.pdb | oxygen_adder > output.pdb
//              oxygen_adder < input.pdb > output.pdb
//              oxygen_adder -i input.pdb -o output.pdb
//

/**
 * A utility function for conversion of strings to everything else.
 */
template<typename T> T parse(std::string const &src) {
    std::istringstream iss(src);
    T rv;
    iss >> rv;
    return rv;
}

/**
 * Trims whitespaces at the beginning and at the end of the given string.
 */
std::string trim(std::string const &src) {
    size_t from = src.find_first_not_of(' ');
    size_t to = src.find_last_not_of(' ');
    if (from == src.npos) {
        return "";
    } else {
        return src.substr(from, to - from + 1);
    }
}

/**
 * A utility which tokenizes an ATOM string into tokens.
 */
void tokenize(std::string const& src, std::vector<std::string> &rv) {
    //ATOM      1  N   TYR A  38     -14.433  26.356  16.716  1.00132.53           N
    //000011111112222233334455556666666666667777777788888888999999000000111111111111
    rv.push_back(trim(src.substr(0, 4)));
    rv.push_back(trim(src.substr(4, 7)));
    rv.push_back(trim(src.substr(11, 5)));
    rv.push_back(trim(src.substr(16, 4)));
    rv.push_back(trim(src.substr(20, 2)));
    rv.push_back(trim(src.substr(22, 4)));
    rv.push_back(trim(src.substr(26, 12)));
    rv.push_back(trim(src.substr(38, 8)));
    rv.push_back(trim(src.substr(46, 8)));
    rv.push_back(trim(src.substr(54, 6)));
    rv.push_back(trim(src.substr(60, 6)));
    rv.push_back(src.substr(66, src.length() - 66));
}

/**
 * An atom which is constructed from parts of ATOM string from PDB.
 */
struct atom {
    std::string type;
    std::string acid;
    std::vector<std::string> extra_tokens;
    int acid_number;
    double x, y, z, t1, t2;

    /**
     * Constructs an atom from a well-formed ATOM string tokenized by spaces.
     */
    atom(std::vector<std::string> const &src) {
        type = src[2];
        acid = src[3];
        acid_number = parse<int>(src[5]);
        x = parse<double>(src[6]);
        y = parse<double>(src[7]);
        z = parse<double>(src[8]);
        t1 = parse<double>(src[9]);
        t2 = parse<double>(src[10]);
        for (unsigned i = 11; i < src.size(); ++i) {
            extra_tokens.push_back(src[i]);
        }
    }
};

/**
 * Prints the given data in the PDB ATOM string format.
 */
void print_atom(unsigned id, std::string const &type, std::string const &acid, unsigned 
acid_number,
                double x, double y, double z, std::ostream &output) {
    fixed(output);
    output << "ATOM";
    output.width(7);
    output << id << "  " << type;
    for (unsigned i = type.length(); i < 4; ++i) output << " ";
    output << acid << " A";
    output.width(4);
    output << acid_number;
    output.width(12);
    output.precision(3);
    output << x;
    output.width(8);
    output.precision(3);
    output << y;
    output.width(8);
    output.precision(3);
    output << z;
    output << "  1.00";
    if (type == "O") {
        output << " 99.99";
    } else {
        output << "  0.00";
    }
    output.width(12);
    output << type[0] << std::endl;
}

double px[] = {1.2345, 6.789};
double py[] = {4.1352, -1.3462};
double pz[] = {0.11414, 5.3552};

/**
 * Does the actual processing of a single model.
 */
int process_atoms(std::vector<atom> &atoms, std::ostream &output) {
    if (check_oxygens) {
        std::vector<atom> filtered;
        for (unsigned i = 0; i < atoms.size(); ++i) {
            std::string const &tp = atoms[i].type;
            if (tp == "N" || tp == "C" || tp == "CA" || tp == "O") {
                filtered.push_back(atoms[i]);
            }
        }
        atoms = filtered;
        std::vector<double> dist_all;
        for (unsigned i = 0; i < atoms.size(); ++i) {
            if (atoms[i].type == "O" && i + 1 < atoms.size()) {
                atom const &me = atoms[i];
                atom const &myC = atoms[i - 1];
                atom const &myCA = atoms[i - 2];
                atom const &nextN = atoms[i + 1];

                double ccx = myC.x - myCA.x;
                double ccy = myC.y - myCA.y;
                double ccz = myC.z - myCA.z;
                double ccl = sqrt(ccx * ccx + ccy * ccy + ccz * ccz);
                ccx /= ccl, ccy /= ccl, ccz /= ccl;

                double ncx = nextN.x - myC.x;
                double ncy = nextN.y - myC.y;
                double ncz = nextN.z - myC.z;
                double ncl = sqrt(ncx * ncx + ncy * ncy + ncz * ncz);
                ncx /= ncl, ncy /= ncl, ncz /= ncl;

                double sm = ccx * ncx + ccy * ncy + ccz * ncz;
                double dx = ccx * sm, dy = ccy * sm, dz = ccz * sm;
                dx += dx - ncx, dy += dy - ncy, dz += dz - ncz;
                double dl = sqrt(dx * dx + dy * dy + dz * dz);
                dx /= dl, dy /= dl, dz /= dl;

                double bond = 1.23;
                double ox = myC.x + dx * bond, oy = myC.y + dy * bond, oz = myC.z + dz * bond;
                dist_all.push_back(sqrt((ox - me.x) * (ox - me.x) + (oy - me.y) * (oy - me.y) + (oz - me.z) * (oz - me.z)));
            }
        }
        std::sort(dist_all.begin(), dist_all.end());
        output << "Differences:" << std::endl;
        for (unsigned i = 0; i < dist_all.size(); ++i) {
            output << dist_all[i] << std::endl;
        }
        double sum = 0;
        for (unsigned i = 0; i < dist_all.size(); ++i) {
            sum += dist_all[i];
        }
        output << "Average: " << sum / dist_all.size() << std::endl;
        output << "Oxygen count: " << dist_all.size() << std::endl;
    } else {
        for (unsigned i = 0; i < atoms.size(); ++i) {
            std::string const &tp = atoms[i].type;
            if (tp != "N" && tp != "C" && tp != "CA") {
                std::cerr << "[oxygen_adder] Unexpected atom role: " << tp << ". This utility works with C, CA, N atoms only." << std::endl;
                return 1;
            }
        }
        unsigned atom_count = 0;
        for (unsigned i = 0; i < atoms.size(); i += 3) {
            atom const &n = atoms[i];
            atom const &ca = atoms[i + 1];
            atom const &c = atoms[i + 2];

            if (n.type  != "N") {
                std::cerr << "[oxygen_adder] Atom " << (i + 1) << ": expected N, found "  << n.type << std::endl;
                return 1;
            }
            if (ca.type != "CA") {
                std::cerr << "[oxygen_adder] Atom " << (i + 2) << ": expected CA, found " << n.type << std::endl;
                return 1;
            }
            if (c.type  != "C") {
                std::cerr << "[oxygen_adder] Atom " << (i + 3) << ": expected C, found "  << n.type << std::endl;
                return 1;
            }

            print_atom(++atom_count,  n.type,  n.acid,  n.acid_number,  n.x,  n.y,  n.z, output);
            print_atom(++atom_count, ca.type, ca.acid, ca.acid_number, ca.x, ca.y, ca.z, output);
            print_atom(++atom_count,  c.type,  c.acid,  c.acid_number,  c.x,  c.y,  c.z, output);

            double ccx = c.x - ca.x;
            double ccy = c.y - ca.y;
            double ccz = c.z - ca.z;
            double ccl = sqrt(ccx * ccx + ccy * ccy + ccz * ccz);
            ccx /= ccl, ccy /= ccl, ccz /= ccl;

            if (i + 3 < atoms.size()) {
                atom const &nn = atoms[i + 3];

                double ncx = nn.x - c.x;
                double ncy = nn.y - c.y;
                double ncz = nn.z - c.z;
                double ncl = sqrt(ncx * ncx + ncy * ncy + ncz * ncz);
                ncx /= ncl, ncy /= ncl, ncz /= ncl;

                double sm = ccx * ncx + ccy * ncy + ccz * ncz;
                double dx = ccx * sm, dy = ccy * sm, dz = ccz * sm;
                dx += dx - ncx, dy += dy - ncy, dz += dz - ncz;
                double dl = sqrt(dx * dx + dy * dy + dz * dz);
                dx /= dl, dy /= dl, dz /= dl;

                double bond = 1.23;
                double ox = c.x + dx * bond, oy = c.y + dy * bond, oz = c.z + dz * bond;
                print_atom(++atom_count, "O", n.acid, n.acid_number, ox, oy, oz, output);
            } else {
                for (unsigned a = 0; a < 2; ++a) {
                    double ppx = px[a], ppy = py[a], ppz = pz[a];
                    double ppl = sqrt(ppx * ppx + ppy * ppy + ppz * ppz);
                    ppx /= ppl, ppy /= ppl, ppz /= ppl;

                    double sm = ccx * ppx + ccy * ppy + ccz * ppz;
                    double devx = ppx - ccx * sm, devy = ppy - ccy * sm, devz = ppz - ccz * sm;
                    double devl = sqrt(devx * devx + devy * devy + devz * devz);
                    if (devl > 0.1) {
                        devx /= devl, devy /= devl, devz /= devl;
                        double xsin = 0.8723;
                        double xcos = 0.4888;
                        double ox = c.x + 1.23 * (ccx * xcos + devx * xsin);
                        double oy = c.y + 1.23 * (ccy * xcos + devy * xsin);
                        double oz = c.z + 1.23 * (ccz * xcos + devz * xsin);
                        print_atom(++atom_count, "O", n.acid, n.acid_number, ox, oy, oz, output);
                        break;
                    }
                }
            }
        }
        return 0;
    }
}

/**
 * Parses the input line by line,
 * collects sequences of ATOM strings and passes it to the 'process_atoms' routine,
 * all other strings are passed to the output as is.
 */
int process(std::istream &input, std::ostream &output) {
    std::string line;
    std::vector<std::string> tokens;
    std::vector<atom> atoms;
    int line_count = 0;
    while (std::getline(input, line)) {
        ++line_count;
        if (line.find("ATOM") == 0) {
            tokens.clear();
            tokenize(line, tokens);
            if (tokens.size() < 11) {
                std::cerr << "[oxygen_adder] Error: in ATOM line " << line_count << " the number of tokens is "
                          << tokens.size() << ", less than 11" << std::endl;
                return 1;
            }
            atoms.push_back(atom(tokens));
        } else {
            if (atoms.size() > 0) {
                int res = process_atoms(atoms, output);
                if (res != 0) {
                    return res;
                }
            }
            atoms.clear();
            output << line << std::endl;
        }
    }
    if (atoms.size() > 0) {
        int res = process_atoms(atoms, output);
        if (res != 0) {
            return res;
        }
    }
}

int main(int argc, char **argv) {
    char *inf = 0, *ouf = 0;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp("-i", argv[i]) && i + 1 < argc) {
            inf = argv[i + 1];
        }
        if (!std::strcmp("-o", argv[i]) && i + 1 < argc) {
            ouf = argv[i + 1];
        }
        if (!std::strcmp("-c", argv[i])) {
            check_oxygens = true;
        }
    }
    std::ifstream *in = inf == 0 ? 0 : new std::ifstream(inf);
    std::ofstream *out = ouf == 0 ? 0 : new std::ofstream(ouf);
    int rv = process(inf == 0 ? std::cin : *in, ouf == 0 ? std::cout : *out);
    if (in) delete in;
    if (out) delete out;
    return rv;
}
