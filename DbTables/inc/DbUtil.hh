#ifndef DbTables_DbUtil_hh
#define DbTables_DbUtil_hh

#include "Offline/DbTables/inc/DbTableCollection.hh"
#include <string>

namespace mu2e {

class DbUtil {
 public:
  static DbTableCollection readFile(std::string const& fn, bool saveCsv = true);
  static void writeFile(std::string const& fn, DbTableCollection const& coll);

  // split a csv string into lines on \n
  static std::vector<std::string> splitCsvLines(std::string const& csv);
  // split a line of csv into columns
  static std::vector<std::string> splitCsv(std::string const& line);
  // clean up a csv row for SQL insert
  static std::string sqlLine(std::string const& line);
  // provide the current local time as a string
  static std::string timeString();
};

}  // namespace mu2e
#endif
