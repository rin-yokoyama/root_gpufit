#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <chrono>
#include "TFile.h"
#include "TTreeReader.h"
#include "TTree.h"
#include "TTreeReaderValue.h"
#include "Interfaces/PulseFitInterface.hpp"

int main(int argc, char **argv)
{
    if (argc < 4)
    {
        std::cout << "Usage: pulseFitTest [input_root_file] [nfit] [nevent]" << std::endl;
        return 1;
    }

    const int kNFit = std::stoi(argv[2]);
    const int kNEvent = std::stoi(argv[3]);

    TFile ofile("scratch/cpu_output.root", "recreate");
    TTree tree("tree", "tree");
    tree.SetDirectory(&ofile);
    std::vector<float> params(5);
    int index = 0;
    int state = 0;
    float chi_square = 0;
    int n_iterations = 0;
    tree.Branch("params", params.data(), "params[5]/F");

    TFile file(argv[1]);
    TTreeReader reader("tree", &file);
    TTreeReaderValue<std::vector<float>> pulse(reader, "pulse");

    PulseFitInterface fitter(kNFit, 504, 0);
    fitter.SetPrepulseRange(20);
    fitter.SetInitialRiseTime(1);
    fitter.SetInitialDecayTime(50);

    auto cpuFit = [&fitter, &tree, &index, &params, &state, &chi_square, &n_iterations]()
    {
        fitter.CallCpufit();
        while (!fitter.ReadResults(index, params, state, chi_square, n_iterations))
        {
            tree.Fill();
        }
    };

    reader.SetEntriesRange(0, kNEvent);
    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    while (reader.Next())
    {
        std::vector<float> p = *pulse;
        if (fitter.AddPulse(p))
        {
            cpuFit();
        }
    }
    cpuFit();
    std::chrono::milliseconds::rep const dt_cpufit = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0).count();
    std::cout << std::fixed << std::setprecision(0);
    std::cout << "Cpufit: nFit = " << kNFit << ", nEvent = " << kNEvent << std::endl;
    std::cout << "time = " << std::setw(13) << static_cast<double>(dt_cpufit) << std::setw(3) << "ms" << std::endl;

    ofile.Write();
    ofile.Close();
    file.Close();

    return 0;
}