//
// Created by André Schnabel on 10.12.17.
//

#pragma once

#include <string>
#include <boost/algorithm/string/join.hpp>
#include "Helpers.h"
#include "json11.hpp"

class EasyCSV {
public:
    explicit EasyCSV(const std::vector<std::string> columnNames, const std::string &_sep = ";");
    explicit EasyCSV(const std::string& headerRow);

    void addRow(const std::vector<double> &cellContents);
    void addRow(const std::vector<std::string> &cellContents);

    std::string toString() const;
    json11::Json toJson() const;

    void persist(const std::string &outFilename) const;

private:
    std::vector<std::string> cols(std::string s) const;

    std::string sep, contentLines, headerLine;
    int ncols;
};