#ifndef SECRECY_DATABASE_H
#define SECRECY_DATABASE_H

#include "relation.h"

#include <string>
#include <vector>
#include <map>
#include <iostream>

// TODO: Move isBColumn to the expression class.
class Database {
    // Relations
    std::vector<std::shared_ptr<Relation>> relations;
    // Foreign Keys
    std::map<std::string, std::string> fks;
    // Name-to-relation mapping
    std::map<std::string, std::shared_ptr<Relation>> index;

    // Metadata
    std::string name;
    std::string versionName;

public:
 
    Database(const std::string &name,
             const std::vector<std::shared_ptr<Relation>> &relations,
             const std::string &version="None");

    ~Database();

    // Add Foreign-Key mappings
    void addFKs(std::map<std::string, std::string>& fks);

    // Returns the database relations
    const std::vector<std::shared_ptr<Relation>>& getRelations() const { return relations; };

    // Returns an index to the database relations
    const std::map<std::string, std::shared_ptr<Relation>>& getIndex() const { return index; };

    /***
     * Retrieves all base relations containing an attribute with name `att_name`.
     * @param att_name - The attribute name.
     * @param relations - The vector of relations to populate.
     */
    void getRelations(const std::string& att_name, std::vector<std::shared_ptr<Relation>>& relations) const;

    /***
     * Retrieves the base relation with name `rel_name`.
     * @param rel_name - The relation name.
     */
    std::shared_ptr<Relation> getRelation(const std::string& rel_name) const;

    bool containsRelation(const std::string& name) const { return index.find(name) != index.end(); };

    unsigned relNum() const { return relations.size(); };

    bool contains(const std::string &att_name, const std::string &rel_name) const;

    // Prints schema to stdout
    void printSchema();
};

#endif //SECRECY_DATABASE_H
