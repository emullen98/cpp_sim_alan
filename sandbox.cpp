#include <iostream>
// Scripts held locally must be put in quotes
#include "cxxopts.hpp"

int main(int argc, char* argv[]) {
    cxxopts::Options options("MyProgram", "A simple program with default arguments.");

    options.add_options()
        ("f,file", "Input file", cxxopts::value<std::string>()->default_value("data.txt"))
        ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
        ("t,time", "Simulation timesteps", cxxopts::value<int>()->default_value("10000"))
        ("s,size", "Simulation size (num of cells)", cxxopts::value<int>()->default_value("10000"))
        ("r,rate", "Strain rate", cxxopts::value<double>()->default_value("0.0"))
        ("d,disorder", "Disorder width", cxxopts::value<double>()->default_value("0.05"))
        ("w,weakening", "Weakening (epsilon)", cxxopts::value<double>()->default_value("0.0"));

    auto result = options.parse(argc, argv);

    bool verbose = result["verbose"].as<bool>();
    int time = result["time"].as<int>();
    int size = result["size"].as<int>();
    double rate = result["rate"].as<double>();
    double disorder = result["disorder"].as<double>();
    double weakening = result["weakening"].as<double>();
    std::string filename = result["file"].as<std::string>();
    std::string out_filename = "stress_s=" + std::to_string(size) + "_r=" + std::to_string(rate) + "_d=" + std::to_string(disorder) + "_w=" + std::to_string(weakening) + ".txt" ;

    std::cout << "Filename in: " << filename << std::endl;
    std::cout << "Verbose: " << verbose << std::endl;
    std::cout << "Simulation time: " << time << std::endl;
    std::cout << "System size: " << size << std::endl;
    std::cout << "Strain rate: " << rate << std::endl;
    std::cout << "Arrest stress distr width: " << disorder << std::endl;
    std::cout << "Weakening (epsilon): " << weakening << std::endl;
    std::cout << "Filename out: " << out_filename << std::endl;

    return 0;
}