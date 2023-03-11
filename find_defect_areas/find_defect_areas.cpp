#include "find_defect_areas.h"
#include <iostream>
#include <fstream>
#include <string>


std::vector<std::vector<float>> preprocessing(const std::vector<std::vector<float>>& thicknesses)
{
    std::vector<std::vector<float>> new_thicknesses = thicknesses;
    for (int i = 0; i < 64; ++i)
    {
        std::vector<float>::iterator last_iter = new_thicknesses[i].end();
        bool flag = false;
        for (auto iter = new_thicknesses[i].begin(); iter != new_thicknesses[i].end(); ++iter)
        {
            if (*iter > 0)
            {
                if (flag)
                {
                    std::fill(last_iter + 1, iter, (*last_iter + *iter) / 2);
                    flag = false;
                }
                last_iter = iter;
            }
            else
            {
                if (last_iter != new_thicknesses[i].end())
                {
                    flag = true;
                }
            }
        }
    }
    return std::move(new_thicknesses);
}

std::vector<float> find_global_medians(const std::vector<std::vector<float>>& thicknesses)
{
    std::vector<float> medians(64);
    for (int i = 0; i < 64; ++i)
    {
        std::vector<float> vec_tmp;
        for (float value : thicknesses[i])
        {
            if (value > 0)
            {
                vec_tmp.push_back(value);
            }
        }
        if (vec_tmp.size() % 2 == 0)
        {
            std::nth_element(vec_tmp.begin(), vec_tmp.begin() + vec_tmp.size() / 2, vec_tmp.end());
            float val_tmp = vec_tmp[vec_tmp.size() / 2];
            std::nth_element(vec_tmp.begin(), vec_tmp.begin() + (vec_tmp.size() - 1) / 2, vec_tmp.end());
            medians[i] = (val_tmp + vec_tmp[(vec_tmp.size() - 1) / 2]) / 2;
        }
        else
        {
            std::nth_element(vec_tmp.begin(), vec_tmp.begin() + vec_tmp.size() / 2, vec_tmp.end());
            medians[i] = vec_tmp[vec_tmp.size() / 2];
        }
    }
    return std::move(medians);
}

float find_local_median(std::vector<float>::iterator start_pos, std::vector<float>::iterator end_pos)
{
    std::vector<float> vec_tmp(end_pos - start_pos);
    std::copy(start_pos, end_pos, vec_tmp.begin());
    std::sort(vec_tmp.begin(), vec_tmp.end());
    if (vec_tmp.size() % 2 == 0)
    {
        return (vec_tmp[vec_tmp.size() / 2] + vec_tmp[(vec_tmp.size() - 1) / 2]) / 2;
    }
    else
    {
        return vec_tmp[vec_tmp.size() / 2];
    }
}

bool check_bounds(std::vector<float>::iterator start_pos, std::vector<float>::iterator end_pos)
{
    for (auto iter = start_pos; iter != end_pos; ++iter)
    {
        if (*iter < 0)
        {
            return false;
        }
    }
    return true;
}

void localizeDefect(const std::vector<std::vector<bool>>& pipe, std::vector<std::vector<bool>>& visitedPositions, std::set<int>& rowIndexSet, std::size_t& minColumnIndex, std::size_t& maxColumnIndex, bool& round_flag, int i, std::size_t j)
{
    visitedPositions[i][j] = true;
    rowIndexSet.insert(i);
    if (j < minColumnIndex)
    {
        minColumnIndex = j;
    }
    else if (j > maxColumnIndex)
    {
        maxColumnIndex = j;
    }
    if (i + 1 < 64 && !visitedPositions[i + 1][j] && pipe[i + 1][j])
    {
        localizeDefect(pipe, visitedPositions, rowIndexSet, minColumnIndex, maxColumnIndex, round_flag, i + 1, j);
    }
    if (i == 63 && !visitedPositions[0][j] && pipe[0][j])
    {
        round_flag = true;
        localizeDefect(pipe, visitedPositions, rowIndexSet, minColumnIndex, maxColumnIndex, round_flag, 0, j);
    }
    if (j + 1 < pipe[i].size() && !visitedPositions[i][j + 1] && pipe[i][j + 1])
    {
        localizeDefect(pipe, visitedPositions, rowIndexSet, minColumnIndex, maxColumnIndex, round_flag, i, j + 1);
    }
    if (j > 0 && !visitedPositions[i][j - 1] && pipe[i][j - 1])
    {
        localizeDefect(pipe, visitedPositions, rowIndexSet, minColumnIndex, maxColumnIndex, round_flag, i, j - 1);
    }
    if (i + 1 < 64 && j + 1 < pipe[i].size() && !visitedPositions[i + 1][j + 1] && pipe[i + 1][j + 1])
    {
        localizeDefect(pipe, visitedPositions, rowIndexSet, minColumnIndex, maxColumnIndex, round_flag, i + 1, j + 1);
    }
    if (i == 63 && j + 1 < pipe[i].size() && !visitedPositions[0][j + 1] && pipe[0][j + 1])
    {
        round_flag = true;
        localizeDefect(pipe, visitedPositions, rowIndexSet, minColumnIndex, maxColumnIndex, round_flag, 0, j + 1);
    }
    if (i + 1 < 64 && j > 0 && !visitedPositions[i + 1][j - 1] && pipe[i + 1][j - 1])
    {
        localizeDefect(pipe, visitedPositions, rowIndexSet, minColumnIndex, maxColumnIndex, round_flag, i + 1, j - 1);
    }
    if (i == 63 && j > 0 && !visitedPositions[0][j - 1] && pipe[0][j - 1])
    {
        round_flag = true;
        localizeDefect(pipe, visitedPositions, rowIndexSet, minColumnIndex, maxColumnIndex, round_flag, 0, j - 1);
    }
    if (i + 2 < 64 && !visitedPositions[i + 2][j] && pipe[i + 2][j])
    {
        localizeDefect(pipe, visitedPositions, rowIndexSet, minColumnIndex, maxColumnIndex, round_flag, i + 2, j);
    }
    if (i + 2 >= 64 && !visitedPositions[(i + 2) % 64][j] && pipe[(i + 2) % 64][j])
    {
        round_flag = true;
        localizeDefect(pipe, visitedPositions, rowIndexSet, minColumnIndex, maxColumnIndex, round_flag, (i + 2) % 64, j);
    }
    if (j + 2 < pipe[i].size() && !visitedPositions[i][j + 2] && pipe[i][j + 2])
    {
        localizeDefect(pipe, visitedPositions, rowIndexSet, minColumnIndex, maxColumnIndex, round_flag, i, j + 2);
    }
    if (j > 1 && !visitedPositions[i][j - 2] && pipe[i][j - 2])
    {
        localizeDefect(pipe, visitedPositions, rowIndexSet, minColumnIndex, maxColumnIndex, round_flag, i, j - 2);
    }
    if (i - 1 >= 0 && !visitedPositions[i - 1][j] && pipe[i - 1][j])
    {
        localizeDefect(pipe, visitedPositions, rowIndexSet, minColumnIndex, maxColumnIndex, round_flag, i - 1, j);
    }
    if (i == 0 && !visitedPositions[63][j] && pipe[63][j])
    {
        round_flag = true;
        localizeDefect(pipe, visitedPositions, rowIndexSet, minColumnIndex, maxColumnIndex, round_flag, 63, j);
    }
    if (i - 1 >= 0 && j + 1 < pipe[i].size() && !visitedPositions[i - 1][j + 1] && pipe[i - 1][j + 1])
    {
        localizeDefect(pipe, visitedPositions, rowIndexSet, minColumnIndex, maxColumnIndex, round_flag, i - 1, j + 1);
    }
    if (i == 0 && j + 1 < pipe[i].size() && !visitedPositions[63][j + 1] && pipe[63][j + 1])
    {
        round_flag = true;
        localizeDefect(pipe, visitedPositions, rowIndexSet, minColumnIndex, maxColumnIndex, round_flag, 63, j + 1);
    }
    if (i - 1 >= 0 && j > 0 && !visitedPositions[i - 1][j - 1] && pipe[i - 1][j - 1])
    {
        localizeDefect(pipe, visitedPositions, rowIndexSet, minColumnIndex, maxColumnIndex, round_flag, i - 1, j - 1);
    }
    if (i == 0 && j > 0 && !visitedPositions[63][j - 1] && pipe[63][j - 1])
    {
        round_flag = true;
        localizeDefect(pipe, visitedPositions, rowIndexSet, minColumnIndex, maxColumnIndex, round_flag, 63, j - 1);
    }
    if (i - 2 >= 0 && !visitedPositions[i - 2][j] && pipe[i - 2][j])
    {
        localizeDefect(pipe, visitedPositions, rowIndexSet, minColumnIndex, maxColumnIndex, round_flag, i - 2, j);
    }
    if (i - 2 < 0 && !visitedPositions[62 + i][j] && pipe[62 + i][j])
    {
        round_flag = true;
        localizeDefect(pipe, visitedPositions, rowIndexSet, minColumnIndex, maxColumnIndex, round_flag, 62 + i, j);
    }
}

bool checkRoundDefect(const std::set<int>& rowIndexSet)
{
    int firstElem;
    if (rowIndexSet.find(0) != rowIndexSet.end())
    {
        firstElem = 0;
    }
    else if (rowIndexSet.find(1) != rowIndexSet.end())
    {
        firstElem = 1;
    }
    else
    {
        return false;
    }
    while (firstElem < 62)
    {
        if (rowIndexSet.find(firstElem + 1) != rowIndexSet.end())
        {
            ++firstElem;
        }
        else if (rowIndexSet.find(firstElem + 2) != rowIndexSet.end())
        {
            firstElem += 2;
        }
        else
        {
            return false;
        }
    }
    if (firstElem == 63)
    {
        return true;
    }
    if (rowIndexSet.find(0) != rowIndexSet.end())
    {
        return true;
    }
    if (rowIndexSet.find(63) != rowIndexSet.end())
    {
        return true;
    }
    return false;
}

std::vector<defect> findDefectZones(const std::vector<std::vector<bool>>& pipe, int min_sensor_coverage, std::size_t start_pos)
{
    std::vector<std::vector<bool>> visitedPositions(64, std::vector<bool>(pipe[0].size(), false));
    std::vector<defect> result;
    for (int i = 0; i < 64; ++i)
    {
        for (std::size_t j = 0; j < pipe[i].size(); ++j)
        {
            if (!visitedPositions[i][j])
            {
                if (!pipe[i][j])
                {
                    visitedPositions[i][j] = true;
                }
                else
                {
                    std::set<int> rowIndexSet;
                    std::size_t minColumnIndex = j, maxColumnIndex = j;
                    bool round_flag = false;
                    localizeDefect(pipe, visitedPositions, rowIndexSet, minColumnIndex, maxColumnIndex, round_flag, i, j);
                    if (checkRoundDefect(rowIndexSet))
                    {
                        defect box;
                        box.start_pos = start_pos + minColumnIndex;
                        box.end_pos = start_pos + maxColumnIndex;
                        box.start_sensor = 0;
                        box.end_sensor = 63;
                        result.push_back(box);
                    }
                    else
                    {
                        if (round_flag)
                        {
                            auto iter = rowIndexSet.begin();
                            int prev = *iter, next = *iter;
                            ++iter;
                            for (; iter != rowIndexSet.end(); ++iter)
                            {
                                next = *iter;
                                if (next - prev > 2)
                                {
                                    break;
                                }
                                prev = next;
                            }
                            if (prev + 65 - next >= min_sensor_coverage)
                            {
                                defect box;
                                box.start_pos = start_pos + minColumnIndex;
                                box.end_pos = start_pos + maxColumnIndex;
                                box.start_sensor = next;
                                box.end_sensor = prev;
                                result.push_back(box);
                            }
                        }
                        else
                        {
                            if (*(rowIndexSet.rbegin()) - *(rowIndexSet.begin()) + 1 >= min_sensor_coverage)
                            {
                                defect box;
                                box.start_pos = start_pos + minColumnIndex;
                                box.end_pos = start_pos + maxColumnIndex;
                                box.start_sensor = *(rowIndexSet.begin());
                                box.end_sensor = *(rowIndexSet.rbegin());
                                result.push_back(box);
                            }
                        }
                    }
                }
            }
        }
    }
    return std::move(result);
}

defects_description find_defect_areas(std::size_t start_pos, const std::vector<std::vector<float>>& thicknesses,
    float max_thickness_deviation_perc, int window_size, int min_sensor_coverage)  // функция вызывается один раз на всей трубе
{
    defects_description defectsDescription;
    defectsDescription.all_defect_positions.resize(64);
    std::vector<std::vector<float>> new_thicknesses = preprocessing(thicknesses);
    std::vector<std::vector<bool>> binary_defects(64, std::vector<bool>(new_thicknesses[0].size(), false));
    std::vector<float> global_medians = find_global_medians(new_thicknesses);
    for (int i = 0; i < 64; ++i)
    {
        std::size_t right_bound = 0;
        for (std::size_t j = 0; j < new_thicknesses[i].size() - window_size + 1; ++j)
        {
            if (check_bounds(new_thicknesses[i].begin() + j, new_thicknesses[i].begin() + j + window_size))
            {
                if (std::abs(find_local_median(new_thicknesses[i].begin() + j, new_thicknesses[i].begin() + j + window_size) - global_medians[i]) > (global_medians[i] * (max_thickness_deviation_perc / 100)))
                {
                    right_bound = j + window_size;
                }
            }
            if (j < right_bound)
            {
                defectsDescription.all_defect_positions[i].push_back(j + start_pos);
                binary_defects[i][j] = true;
            }
        }
        for (std::size_t j = new_thicknesses[i].size() - window_size + 1; j < right_bound; ++j)
        {
            defectsDescription.all_defect_positions[i].push_back(j + start_pos);
            binary_defects[i][j] = true;
        }
    }
    defectsDescription.defects = findDefectZones(binary_defects, min_sensor_coverage, start_pos);
    return std::move(defectsDescription);
}

int main()
{
    // real pipe test

    std::size_t start_pos = 7802;
    std::vector<std::vector<float>> thicknesses(64);
    std::ifstream fin;

    float a;
    for (int i = 0; i < 64; ++i)
    {
        fin.open("data/" + std::to_string(i) + ".txt");
        std::vector<float> vec;
        while (fin)
        {
            fin >> a;
            vec.push_back(a);
        }
        vec.pop_back();
        thicknesses[i] = vec;
        fin.close();
    }

    defects_description defectsDescription = find_defect_areas(start_pos, thicknesses, 30.0, 4, 3);

    for (int i = 0; i < 64; ++i)
    {
        if (defectsDescription.all_defect_positions[i].size() > 0)
        {
            std::cout << "Defects at sensor: " + std::to_string(i) << std::endl;
            for (std::size_t pos : defectsDescription.all_defect_positions[i])
            {
                std::cout << pos << ' ';
            }
            std::cout << std::endl;
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
    std::cout << "Isolated defects found: " + std::to_string(defectsDescription.defects.size()) << std::endl;
    for (const defect& def : defectsDescription.defects)
    {
        std::cout << "Start pos: " + std::to_string(def.start_pos) << std::endl;
        std::cout << "End pos: " + std::to_string(def.end_pos) << std::endl;
        std::cout << "Start sensor: " + std::to_string(def.start_sensor) << std::endl;
        std::cout << "End sensor: " + std::to_string(def.end_sensor) << std::endl;
        std::cout << std::endl;
    }

    // synthetic round test
    std::cout << "Synthetic round test" << std::endl;

    thicknesses.clear();
    thicknesses.resize(64);
    start_pos = 0;
    for (int i = 0; i < 64; ++i)
    {
        std::vector<float> vec_tmp(200);
        for (int j = 0; j < 200; ++j)
        {
            if ((i < 2 || i == 63) && j > 100 && j < 110)
            {
                vec_tmp[j] = 2.0;
            }
            else
            {
                vec_tmp[j] = 5.0;
            }
        }
        thicknesses[i] = vec_tmp;
    }

    defectsDescription = find_defect_areas(start_pos, thicknesses, 30.0, 4, 3);

    for (int i = 0; i < 64; ++i)
    {
        if (defectsDescription.all_defect_positions[i].size() > 0)
        {
            std::cout << "Defects at sensor: " + std::to_string(i) << std::endl;
            for (std::size_t pos : defectsDescription.all_defect_positions[i])
            {
                std::cout << pos << ' ';
            }
            std::cout << std::endl;
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
    std::cout << "Isolated defects found: " + std::to_string(defectsDescription.defects.size()) << std::endl;
    for (const defect& def : defectsDescription.defects)
    {
        std::cout << "Start pos: " + std::to_string(def.start_pos) << std::endl;
        std::cout << "End pos: " + std::to_string(def.end_pos) << std::endl;
        std::cout << "Start sensor: " + std::to_string(def.start_sensor) << std::endl;
        std::cout << "End sensor: " + std::to_string(def.end_sensor) << std::endl;
        std::cout << std::endl;
    }

    // synthetic round test 2

    std::cout << "Synthetic round test 2" << std::endl;

    thicknesses.clear();
    thicknesses.resize(64);
    start_pos = 0;
    for (int i = 0; i < 64; ++i)
    {
        std::vector<float> vec_tmp(200);
        for (int j = 0; j < 200; ++j)
        {
            vec_tmp[j] = 5.0;
        }
        thicknesses[i] = vec_tmp;
    }

    for (int i = 0; i < 4; ++i)
    {
        for (int j = 100; j < 105; ++j)
        {
            thicknesses[i][j] = 2.0;
        }
    }
    for (int j = 105; j < 115; ++j)
    {
        thicknesses[3][j] = 2.0;
    }
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 115; j < 120; ++j)
        {
            thicknesses[i][j] = 2.0;
        }
    }
    for (int i = 62; i < 64; ++i)
    {
        for (int j = 100; j < 105; ++j)
        {
            thicknesses[i][j] = 2.0;
        }
        for (int j = 115; j < 120; ++j)
        {
            thicknesses[i][j] = 2.0;
        }
    }

    defectsDescription = find_defect_areas(start_pos, thicknesses, 30.0, 4, 3);

    for (int i = 0; i < 64; ++i)
    {
        if (defectsDescription.all_defect_positions[i].size() > 0)
        {
            std::cout << "Defects at sensor: " + std::to_string(i) << std::endl;
            for (std::size_t pos : defectsDescription.all_defect_positions[i])
            {
                std::cout << pos << ' ';
            }
            std::cout << std::endl;
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
    std::cout << "Isolated defects found: " + std::to_string(defectsDescription.defects.size()) << std::endl;
    for (const defect& def : defectsDescription.defects)
    {
        std::cout << "Start pos: " + std::to_string(def.start_pos) << std::endl;
        std::cout << "End pos: " + std::to_string(def.end_pos) << std::endl;
        std::cout << "Start sensor: " + std::to_string(def.start_sensor) << std::endl;
        std::cout << "End sensor: " + std::to_string(def.end_sensor) << std::endl;
        std::cout << std::endl;
    }

    // synthetic round test 3

    std::cout << "Synthetic round test 3" << std::endl;

    thicknesses.clear();
    thicknesses.resize(64);
    start_pos = 0;
    for (int i = 0; i < 64; ++i)
    {
        std::vector<float> vec_tmp(200);
        for (int j = 0; j < 200; ++j)
        {
            vec_tmp[j] = 5.0;
        }
        thicknesses[i] = vec_tmp;
    }

    for (int i = 0; i < 64; ++i)
    {
        if (i % 2 == 0)
        {
            for (int j = 100; j < 105; ++j)
            {
                thicknesses[i][j] = 2.0;
            }
        }
    }

    defectsDescription = find_defect_areas(start_pos, thicknesses, 30.0, 4, 64);

    for (int i = 0; i < 64; ++i)
    {
        if (defectsDescription.all_defect_positions[i].size() > 0)
        {
            std::cout << "Defects at sensor: " + std::to_string(i) << std::endl;
            for (std::size_t pos : defectsDescription.all_defect_positions[i])
            {
                std::cout << pos << ' ';
            }
            std::cout << std::endl;
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
    std::cout << "Isolated defects found: " + std::to_string(defectsDescription.defects.size()) << std::endl;
    for (const defect& def : defectsDescription.defects)
    {
        std::cout << "Start pos: " + std::to_string(def.start_pos) << std::endl;
        std::cout << "End pos: " + std::to_string(def.end_pos) << std::endl;
        std::cout << "Start sensor: " + std::to_string(def.start_sensor) << std::endl;
        std::cout << "End sensor: " + std::to_string(def.end_sensor) << std::endl;
        std::cout << std::endl;
    }

    // synthetic test 4

    std::cout << "Synthetic test 4" << std::endl;

    thicknesses.clear();
    thicknesses.resize(64);
    start_pos = 0;
    for (int i = 0; i < 64; ++i)
    {
        std::vector<float> vec_tmp(200);
        for (int j = 0; j < 200; ++j)
        {
            vec_tmp[j] = 5.0;
        }
        thicknesses[i] = vec_tmp;
    }

    thicknesses[4][25] = 2.0;
    thicknesses[4][26] = 1.0;
    thicknesses[4][27] = 2.5;
    thicknesses[4][28] = 2.0;
    thicknesses[5][25] = 2.0;
    thicknesses[5][26] = 1.0;
    thicknesses[5][27] = 2.5;
    thicknesses[5][28] = 2.0;

    thicknesses[14][25] = 2.0;
    thicknesses[14][26] = 1.0;
    thicknesses[14][27] = 2.5;
    thicknesses[14][28] = 2.0;
    thicknesses[15][25] = 2.0;
    thicknesses[15][26] = 1.0;
    thicknesses[15][27] = 2.5;
    thicknesses[15][28] = 2.0;
    thicknesses[17][25] = 2.0;
    thicknesses[17][26] = 1.0;
    thicknesses[17][27] = 2.5;
    thicknesses[17][28] = 2.0;
    thicknesses[17][30] = 2.0;
    thicknesses[17][31] = 1.0;
    thicknesses[17][32] = 2.5;
    thicknesses[17][33] = 2.0;

    thicknesses[2][25] = 2.0;
    thicknesses[2][26] = 1.0;
    thicknesses[2][27] = 2.5;
    thicknesses[2][28] = 2.0;
    thicknesses[0][25] = 2.0;
    thicknesses[0][26] = 1.0;
    thicknesses[0][27] = 2.5;
    thicknesses[0][28] = 2.0;
    thicknesses[62][25] = 2.0;
    thicknesses[62][26] = 1.0;
    thicknesses[62][27] = 2.5;
    thicknesses[62][28] = 2.0;
    thicknesses[61][28] = 2.0;
    thicknesses[61][29] = 1.0;
    thicknesses[61][30] = 2.5;
    thicknesses[61][31] = 2.0;

    defectsDescription = find_defect_areas(start_pos, thicknesses, 30.0, 4, 2);

    for (int i = 0; i < 64; ++i)
    {
        if (defectsDescription.all_defect_positions[i].size() > 0)
        {
            std::cout << "Defects at sensor: " + std::to_string(i) << std::endl;
            for (std::size_t pos : defectsDescription.all_defect_positions[i])
            {
                std::cout << pos << ' ';
            }
            std::cout << std::endl;
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
    std::cout << "Isolated defects found: " + std::to_string(defectsDescription.defects.size()) << std::endl;
    for (const defect& def : defectsDescription.defects)
    {
        std::cout << "Start pos: " + std::to_string(def.start_pos) << std::endl;
        std::cout << "End pos: " + std::to_string(def.end_pos) << std::endl;
        std::cout << "Start sensor: " + std::to_string(def.start_sensor) << std::endl;
        std::cout << "End sensor: " + std::to_string(def.end_sensor) << std::endl;
        std::cout << std::endl;
    }

    // synthetic round test 5

    std::cout << "Synthetic round test 5" << std::endl;

    thicknesses.clear();
    thicknesses.resize(64);
    start_pos = 0;
    for (int i = 0; i < 64; ++i)
    {
        std::vector<float> vec_tmp(200);
        for (int j = 0; j < 200; ++j)
        {
            vec_tmp[j] = 5.0;
        }
        thicknesses[i] = vec_tmp;
    }

    for (int i = 60; i < 64; ++i)
    {
        for (int j = 100; j < 105; ++j)
        {
            thicknesses[i][j] = 2.0;
        }
    }
    for (int j = 105; j < 115; ++j)
    {
        thicknesses[60][j] = 2.0;
    }
    for (int i = 60; i < 64; ++i)
    {
        for (int j = 115; j < 120; ++j)
        {
            thicknesses[i][j] = 2.0;
        }
    }
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 100; j < 105; ++j)
        {
            thicknesses[i][j] = 2.0;
        }
        for (int j = 115; j < 125; ++j)
        {
            thicknesses[i][j] = 2.0;
        }
    }

    defectsDescription = find_defect_areas(start_pos, thicknesses, 30.0, 4, 3);

    for (int i = 0; i < 64; ++i)
    {
        if (defectsDescription.all_defect_positions[i].size() > 0)
        {
            std::cout << "Defects at sensor: " + std::to_string(i) << std::endl;
            for (std::size_t pos : defectsDescription.all_defect_positions[i])
            {
                std::cout << pos << ' ';
            }
            std::cout << std::endl;
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
    std::cout << "Isolated defects found: " + std::to_string(defectsDescription.defects.size()) << std::endl;
    for (const defect& def : defectsDescription.defects)
    {
        std::cout << "Start pos: " + std::to_string(def.start_pos) << std::endl;
        std::cout << "End pos: " + std::to_string(def.end_pos) << std::endl;
        std::cout << "Start sensor: " + std::to_string(def.start_sensor) << std::endl;
        std::cout << "End sensor: " + std::to_string(def.end_sensor) << std::endl;
        std::cout << std::endl;
    }
    return 0;
}
