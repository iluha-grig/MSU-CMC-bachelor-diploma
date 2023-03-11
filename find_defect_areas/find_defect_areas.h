#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <set>


struct defect
{
    std::size_t start_pos;
    std::size_t end_pos;
    int start_sensor;
    int end_sensor;
};

struct defects_description
{
    std::vector<std::vector<std::size_t>> all_defect_positions;
    std::vector<defect> defects;
};

defects_description find_defect_areas(std::size_t start_pos, const std::vector<std::vector<float>>& thicknesses,
    float max_thickness_deviation_perc = 30.0, int window_size = 4,
    int min_sensor_coverage = 3);
