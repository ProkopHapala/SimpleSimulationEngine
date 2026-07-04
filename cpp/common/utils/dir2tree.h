#ifndef dir2tree_h
#define dir2tree_h

/// @file dir2tree.h
/// @brief Populate a TreeView from the filesystem. Depends on GUI.h (TreeViewTree).

#include <dirent.h>
#include <sys/stat.h>
#include <algorithm>
#include <string>
#include <vector>
#include <functional>

inline void dir2tree(TreeViewTree& node, const char* name, const std::string& prefix = "",
                     bool bShowAll = false,
                     std::function<bool(const std::string&)> fileFilter = nullptr) {
    node.content.caption = name;
    std::string path = (prefix.length() == 0) ? name : (prefix + "/" + name);
    node.content.path = path;
    DIR* dir = opendir(path.c_str());
    if (dir) {
        node.content.bDir = true;
        std::vector<std::pair<std::string, bool>> entries;
        struct dirent* ent;
        while ((ent = readdir(dir)) != nullptr) {
            if (ent->d_name[0] == '.') continue;
            std::string childName = ent->d_name;
            std::string childPath = path + "/" + childName;
            struct stat st;
            if (stat(childPath.c_str(), &st) != 0) continue;
            bool isDir = S_ISDIR(st.st_mode);
            if (!isDir && !bShowAll && fileFilter && !fileFilter(childName)) continue;
            entries.push_back({childName, isDir});
        }
        closedir(dir);
        std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
            if (a.second != b.second) return a.second > b.second;
            return a.first < b.first;
        });
        for (auto& [childName, isDir] : entries) {
            TreeViewTree* child = new TreeViewTree();
            child->parrent = &node;
            node.branches.push_back(child);
            dir2tree(*child, childName.c_str(), path, bShowAll, fileFilter);
        }
    } else {
        node.content.bDir = false;
    }
}

#endif
