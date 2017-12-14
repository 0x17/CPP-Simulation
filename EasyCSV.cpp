//
// Created by André Schnabel on 10.12.17.
//

#include <boost/algorithm/string.hpp>
#include "EasyCSV.h"

using namespace std;

EasyCSV::EasyCSV(const vector<string> columnNames, const string &_sep) : sep(_sep), contentLines(""), headerLine(boost::algorithm::join(columnNames, sep)), ncols(columnNames.size()) {}

EasyCSV::EasyCSV(const string &headerRow) : contentLines(""), headerLine(headerRow), ncols(std::count(headerRow.begin(), headerRow.end(), ';')+1), sep(";") {}

void EasyCSV::addRow(const vector<double> &cellContents) {
    addRow(Helpers::constructVector<string>(cellContents.size(), [&cellContents](int i) { return to_string(cellContents[i]); }));
}

void EasyCSV::addRow(const vector<string> &cellContents) {
    contentLines += boost::algorithm::join(cellContents, sep)  + "\n";
}

string EasyCSV::toString() const {
    return headerLine + "\n" + contentLines;
}

void EasyCSV::persist(const string &outFilename) const {
    Helpers::spit(toString(), outFilename + ".csv");
    Helpers::spit(toJson().dump(), outFilename + ".json");
}

vector<string> EasyCSV::cols(string s) const {
    vector<string> res(ncols);
    boost::algorithm::split(res, s, boost::is_any_of(";"));
    return res;
}

json11::Json EasyCSV::toJson() const {
    vector<string> lines;
    boost::split(lines,contentLines,boost::is_any_of("\n"));
	if (lines[lines.size() - 1].empty())
		lines.pop_back();

    vector<json11::Json> rows = Helpers::constructVector<json11::Json>(lines.size(), [this, &lines](int rowIndex) {
        return cols(lines[rowIndex]);
    });

    return json11::Json::object {
            { "header", cols(headerLine) },
            { "rows", rows }
    };
}
