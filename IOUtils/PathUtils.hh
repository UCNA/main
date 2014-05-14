#ifndef PATHUTILS_HH
#define PATHUTILS_HH

#include <string>
#include <vector>

/// check if file exists
bool fileExists(std::string f);
/// check if directory exists
bool dirExists(std::string d);
/// make sure the specified path exists (if not, create it); optionally, exclude last item on path (filename)
void makePath(std::string p, bool forFile = false);
/// list directory contents
std::vector<std::string> listdir(const std::string& dir, bool includeHidden = false);
/// get time since last file modification (s)
double fileAge(const std::string& fname);
/// get environment variable, with default or fail if missing
std::string getEnvSafe(const std::string& v, const std::string& dflt = "FAIL_IF_MISSING");

#endif
