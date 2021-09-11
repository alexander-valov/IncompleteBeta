#include <fstream>
#include <sstream>
#include <iostream>

#include "doctest.h"
#include "JustMath/Beta.hpp"

TEST_SUITE("CompareTable") {

TEST_CASE("CompareTable.CompareTableDouble") {
    std::ifstream table_stream("incbeta_table.csv");

    std::string line;
    int index = 1;
    while(std::getline(table_stream, line)) {
        std::istringstream s(line);
        std::string field;

        /* read x */
        std::getline(s, field, ',');
        double x = std::stod(field);

        /* read a */
        std::getline(s, field, ',');
        double a = std::stod(field);

        /* read b */
        std::getline(s, field, ',');
        double b = std::stod(field);

        /* read res_exact */
        std::getline(s, field, ',');
        double res_exact = std::stod(field);

        CHECK(JustMath::incbeta(x, a, b) == doctest::Approx(res_exact));
        index++;
    }

}

TEST_CASE("CompareTable.CompareTableFloat") {
    std::ifstream table_stream("incbeta_table.csv");

    std::string line;
    int index = 1;
    while(std::getline(table_stream, line)) {
        std::istringstream s(line);
        std::string field;

        /* read x */
        std::getline(s, field, ',');
        float x = std::stof(field);

        /* read a */
        std::getline(s, field, ',');
        float a = std::stof(field);

        /* read b */
        std::getline(s, field, ',');
        float b = std::stof(field);

        /* read res_exact */
        std::getline(s, field, ',');
        float res_exact = std::stod(field);

        CHECK(JustMath::incbeta(x, a, b) == doctest::Approx(res_exact).scale(5));
        index++;
    }

}

}