#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <string>
#include <map>
#include <iomanip>
#include <chrono>

struct BPParameter {
    double L;
    double H;
    BPParameter(double l, double h) : L(l), H(h) {}
};

struct Job {
    int arrival_time;
    int job_size;
};

class JobGenerator {
private:
    std::mt19937 rng;
    const double ALPHA = 1.1;

    // Helper function to generate Pareto distributed random numbers
    int generateParetoSample(double xmin, double xmax) {
        std::uniform_real_distribution<double> u(0.0, 1.0);
        double alpha = ALPHA;
        double u_sample = u(rng);
        
        // Inverse transform sampling for bounded Pareto
        double x = xmin * std::pow(
            (1 - u_sample + u_sample * std::pow(xmin/xmax, alpha)),
            -1/alpha
        );
        return std::ceil(x);
    }

    // Helper function to generate exponential random numbers
    int generateExponentialSample(double mean) {
        std::exponential_distribution<double> exp(1.0/mean);
        return std::max(1, static_cast<int>(std::round(exp(rng))));
    }

public:
    JobGenerator() : rng(std::random_device{}()) {}

    std::vector<Job> generateJobs(int num_jobs, int avg_inter_arrival_time, 
                                 double xmin, double xmax) {
        std::vector<Job> jobs;
        std::vector<int> job_sizes;

        // Generate job sizes
        while (job_sizes.size() < num_jobs) {
            int size = generateParetoSample(xmin, xmax);
            if (size >= xmin && size <= xmax) {
                job_sizes.push_back(size);
            }
        }

        // Generate arrival times
        int current_time = 0;
        std::vector<int> arrival_times;
        for (int i = 0; i < num_jobs; i++) {
            int inter_arrival = generateExponentialSample(avg_inter_arrival_time);
            current_time += inter_arrival;
            arrival_times.push_back(current_time);
        }

        // Combine into jobs
        for (int i = 0; i < num_jobs; i++) {
            jobs.push_back({arrival_times[i], job_sizes[i]});
        }

        return jobs;
    }

    void saveToCSV(const std::string& filename, const std::vector<Job>& jobs) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        // Write header
        file << "arrival_time,job_size\n";

        // Write data
        for (const auto& job : jobs) {
            file << job.arrival_time << "," << job.job_size << "\n";
        }

        file.close();
    }
};

int main() {
    // Define parameters
    std::vector<int> inter_arrival_time;
    for (int i = 20; i <= 40; i += 2) {
        inter_arrival_time.push_back(i);
    }

    std::vector<BPParameter> bp_parameter = {
        {16.772, std::pow(2, 6)},
        {7.918, std::pow(2, 9)},
        {5.649, std::pow(2, 12)},
        {4.639, std::pow(2, 15)},
        {4.073, std::pow(2, 18)}
    };

    JobGenerator generator;
    int num_jobs = 100000;

    // Progress tracking
    int total_combinations = inter_arrival_time.size() * bp_parameter.size();
    int current_combination = 0;

    for (int avg_inter_arrival : inter_arrival_time) {
        for (const auto& bp : bp_parameter) {
            // Progress output
            current_combination++;
            std::cout << "Processing combination " << current_combination 
                      << " of " << total_combinations << " ("
                      << (current_combination * 100.0 / total_combinations) 
                      << "%)\n";

            // Generate jobs
            auto jobs = generator.generateJobs(num_jobs, avg_inter_arrival, 
                                            bp.L, bp.H);

            // Create filename
            std::string filename = "data/(" + 
                                 std::to_string(avg_inter_arrival) + ", " + 
                                 std::to_string(static_cast<int>(bp.L)) + 
                                 ").csv";

            // Save to file
            generator.saveToCSV(filename, jobs);
        }
    }

    return 0;
}