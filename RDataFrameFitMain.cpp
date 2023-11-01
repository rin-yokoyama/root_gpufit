#include <iostream>
#include <string>
#include <numeric>
#include <vector>
#include <iomanip>
#include <chrono>
#include "TROOT.h"
#include "TGraph.h"
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#include "ROOT/RDataFrame.hxx"
#include "TUUID.h"
#include "TMinuit.h"

// Number of workers
const UInt_t kNWorkers = 64U;

// Function to return average value of the vector
template <typename T>
double getAverage(std::vector<T> const &v)
{
    if (v.empty())
    {
        return 0;
    }
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

// Main function
int main(int argc, char **argv)
{
    // Enable multithreading
    // ROOT::EnableImplicitMT(kNWorkers);
    ROOT::EnableImplicitMT();
    ROOT::EnableThreadSafety();
    // Create RDataFrame from a tree, "tree" in the argv[1] file.
    ROOT::RDataFrame d("tree", argv[1]);

    // const int kNEvent = std::stoi(argv[2]);
    auto kNEvent = 2396090;

    // Function to define initial parameters for fitting
    const auto initialParam = [](const std::vector<float> &input)
    {
        std::vector<float> iparams(5);
        const int prePulseRange = 20;
        const float riseTime = 1.0;
        const float decayTime = 57.0;
        const float peakTime = 30.0;

        // Calculate baseline
        const std::vector<float> pre_pulse = {input.begin(), input.begin() + prePulseRange};
        const auto baseline = getAverage(pre_pulse);

        // Calculate amplitude
        const std::vector<float>::const_iterator min_itr = std::min_element(input.begin(), input.end());
        const float amplitude = *min_itr - baseline;

        iparams[0] = amplitude;
        iparams[1] = peakTime;
        iparams[2] = riseTime;
        iparams[3] = decayTime;
        iparams[4] = baseline;
        return iparams;
    };

    // Function to define sample id
    const auto sampleId = [](const std::vector<float> &input)
    {
        std::vector<float> id(input.size());
        std::iota(id.begin(), id.end(), 1);
        return id;
    };

    // Function to fit
    const auto fit = [](const std::vector<float> &input, const std::vector<float> &init_params, const std::vector<float> &id)
    {
        std::vector<float> params(11, 0);
        if (input.empty()||input.size()!=id.size()||init_params.size()!=5)
            return params;
        TUUID uuid;
        auto f1 = new TF1(uuid.AsString(), "0.5*[0]*(1+TMath::Erf((x-[1])/[2]))*exp(-(x-[1])/[3])+[4]", 0, 504);
        f1->SetParameters(init_params[0], init_params[1], init_params[2], init_params[3], init_params[4]);
        TGraph gr(input.size(), id.data(), input.data());
        gr.Fit(f1, "Q", "", 0, 500);
        params[0] = f1->GetParameter(0);
        params[1] = f1->GetParameter(1);
        params[2] = f1->GetParameter(2);
        params[3] = f1->GetParameter(3);
        params[4] = f1->GetParameter(4);
        params[5] = f1->GetParError(0);
        params[6] = f1->GetParError(1);
        params[7] = f1->GetParError(2);
        params[8] = f1->GetParError(3);
        params[9] = f1->GetParError(4);
        params[10] = f1->GetChisquare();
        delete f1;
        return params;
    };

    // auto filtered = d.Range(0,kNEvent);

    std::cout << "Scan started..." << std::endl;
    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    auto output = d.Define("sampleId", sampleId, {"pulse"})
                      .Define("initialParams", initialParam, {"pulse"})
                      .Define("fitParams", fit, {"pulse", "initialParams", "sampleId"});

    output.Snapshot("output", "scratch/rdframe_output.root", {"fitParams"});

    std::chrono::milliseconds::rep const dt_rdframefit = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0).count();
    std::cout << std::fixed << std::setprecision(0);
    std::cout << "RDataFrameFit: nEvent = " << kNEvent << std::endl;
    std::cout << "time = " << std::setw(13) << static_cast<double>(dt_rdframefit) << std::setw(3) << "ms" << std::endl;

    return 0;
}