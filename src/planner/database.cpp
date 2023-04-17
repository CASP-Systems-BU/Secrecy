#include "../../include/planner/database.h"


Database::Database(const std::string &name, 
                    const std::vector<std::shared_ptr<Relation>> &relations,
                    const std::string &version) :
        name(name), relations(relations), versionName(version) {
            // Create relation index
            for (auto r : this->relations) {
                auto ret = this->index.insert(std::pair<std::string,
                                                          std::shared_ptr<Relation>> (r->getName(), r));
                if (ret.second == false) {
                    std::cout << "[ERROR]: Duplicate relation found." << std::endl;
                    exit(-1);
                }
            }
}

Database::~Database() {}

void Database::addFKs(std::map<std::string, std::string> &fks) {
    this->fks = fks;
}

// Retrieves the base relations containing an attribute with name `att_name` (if any)
void Database::getRelations(const std::string& att_name, std::vector<std::shared_ptr<Relation>>& relations) const {
    for (auto r : this->relations) {
        auto cols = r->getCols(true);
        for (auto c : cols) {
            if (c == att_name)
                relations.push_back(r);
        }
    }
}

bool Database::contains(const std::string &att_name, const std::string &rel_name) const {
    for (auto r : this->relations) {
        if (r->getName() == rel_name) {
            auto cols = r->getCols(true);
            for (auto c : cols) {
                if (c == att_name)
                    return true;
            }
        }
    }
    return false;
}

// Returns the base relation with the given name (if any)
std::shared_ptr<Relation> Database::getRelation(const std::string& rel_name) const {
    for (auto r : this->relations) {
        if (r->getName() == rel_name)
            return r;
    }
    return nullptr;
}

void Database::printSchema() {
    if (get_rank()==0) {
        std::cout << "DB Schema: " << std::endl;
        for (int i=0; i<relations.size(); i++) {
            std::cout << relations[i]->getName() << "(";
            auto cols = relations[i]->getCols(true);
            for (int j=0; j<cols.size()-1; j++)
                std::cout << cols[j] << ",";
            std::cout << cols[cols.size()-1] << ")" << std::endl;
        }
    }
}
