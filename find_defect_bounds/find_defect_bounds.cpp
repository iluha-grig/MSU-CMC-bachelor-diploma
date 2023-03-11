#include "find_defect_bounds.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>


int num_resonance(float* spectrum)
{
    float lower_bound = 0.2f, upper_bound = 0.8f;
    int counter = 0;
    bool resonance_flag = false;
    if (spectrum[0] > upper_bound)
    {
        resonance_flag = true;
        ++counter;
    }
    for (int i = 1; i < 249; ++i)
    {
        if (resonance_flag)
        {
            if (spectrum[i] < lower_bound)
            {
                resonance_flag = false;
            }
        }
        else
        {
            if (spectrum[i] > upper_bound)
            {
                resonance_flag = true;
                ++counter;
            }
        }
    }
    return counter;
}

std::vector<float> find_defect_bounds(const std::vector<std::vector<float*>>& spectrums)
{
    std::vector<float> probabilities(64);
    for (int i = 0; i < 64; ++i)
    {
        int defect_spectrum_counter = 0;
        for (std::size_t j = 0; j < spectrums[i].size(); ++j)
        {
            if (num_resonance(spectrums[i][j]) > 1)
            {
                ++defect_spectrum_counter;
            }
        }
        probabilities[i] = static_cast<float>(defect_spectrum_counter) / spectrums[i].size();
    }
    return std::move(probabilities);
}

int main()
{
    // real pipe test

    std::vector<std::vector<float*>> spectrums(64);
    std::ifstream fin;
    for (int i = 0; i < 64; ++i)
    {
        fin.open("data/" + std::to_string(i) + ".txt");
        std::vector<float*> vec_tmp;

        std::vector<std::string> vec;
        while (fin)
        {
            std::string str;
            std::getline(fin, str);
            vec.push_back(str);
        }
        vec.pop_back();
        fin.close();

        for (const std::string& str : vec)
        {
            std::stringstream ss(str);
            float* ptr = new float[249];
            for (int i = 0; i < 249; ++i)
            {
                ss >> ptr[i];
            }
            vec_tmp.push_back(ptr);
        }

        spectrums[i] = vec_tmp;
    }
    
    std::vector<float> probabilities = find_defect_bounds(spectrums);
    for (int i = 0; i < 64; ++i)
    {
        std::cout << "Probability that position at " + std::to_string(i) + " sensor is defect bound = " + std::to_string(probabilities[i]) << std::endl;
    }

    // test num_resonance

    fin.open("data/spectrum_2resonances.txt");
    float* arr = new float[249];
    int i = 0;
    float a;
    while (fin)
    {
        fin >> a;
        if (i < 249)
        {
            arr[i++] = a;
        }
    }
    fin.close();
    std::cout << num_resonance(arr) << std::endl;  // should be 2

    fin.open("data/spectrum_1resonance.txt");
    i = 0;
    while (fin)
    {
        fin >> a;
        if (i < 249)
        {
            arr[i++] = a;
        }
    }
    fin.close();
    std::cout << num_resonance(arr) << std::endl;  // should be 1
    return 0;
}
