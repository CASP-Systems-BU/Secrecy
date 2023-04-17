#include "../../include/planner/relation.h"

Relation::Relation(const std::string &name,
                   const std::vector<std::string> &cols,
                   const int &rows) :
        name(name), materialized(false), aSelectionUsed(false), bSelectionUsed(false) {

    // TODO: use the IDs
    shareTable.numRows = rows;
    shareTable.numCols = 0;
    bShares.resize(cols.size(), false);
    for (int i=0; i<cols.size(); i++) {
        if (isBColumnName(cols[i])) {
            bShares[i] = true;
        }
    }
    addColumns(cols);
}

Relation::Relation(const std::string &name,
                   const std::vector<std::string> &cols,
                   const int &rows,
                   std::vector<bool> &sTypes) :
        name(name), materialized(false), aSelectionUsed(false), bSelectionUsed(false) {

    // TODO: use the IDs
    shareTable.numRows = rows;
    shareTable.numCols = 0;
    colNames = cols;
    bShares = sTypes;
    assert(colNames.size() == bShares.size());
    addColumns(cols);
}

Relation::Relation() {}
Relation::~Relation() {}

void Relation::freeRelation() {
    if (materialized) {
        freeShareTable<ShareTable>(shareTable);
    }

    for (auto it = schema.begin(); it != schema.end(); it++) {
        free(it->second);
    }
}

void Relation::addColumns(const std::vector<std::string> &cols) {
    for (const std::string &col : cols) {
        addColumns(col);
    }
}

void Relation::addColumns(const std::string &col) {
    int index = shareTable.numCols;
    auto c = DerivedColumn<Share, ShareTable>(&shareTable, index);
    auto idx = std::find(colNames.begin(), colNames.end(), col) - colNames.begin();
    assert(idx < bShares.size());
    if (isBColumnName(col) || bShares[idx]) {
        schema[col] = new BColumn<Share, ShareTable>(c);
    } else {
        schema[col] = new AColumn<Share, ShareTable>(c);
    }
    shareTable.numCols++;
}


void Relation::allocateTable() {
    if (!materialized) {
        allocateContent<Share, ShareTable>(shareTable);

        for (auto it = schema.begin(); it != schema.end(); it++) {
            auto other_ = reinterpret_cast<DerivedColumn<Share, ShareTable> *>(it->second);
            other_->allocateRowShares(&shareTable);
        }

        materialized = true;
    }
}

void Relation::insertTable(const ShareTable &other) {
    allocateTable();

    int max_rows = shareTable.numRows < other.numRows ? shareTable.numRows : other.numRows;
    int max_cols = shareTable.numCols < other.numCols ? shareTable.numCols : other.numCols;

    if (shareTable.content != nullptr && other.content != nullptr) {
        for (int i = 0; i < max_rows; ++i) {
            for (int j = 0; j < CR_P.PARTY_REPLICATION * max_cols; ++j) {
                shareTable.content[i][j] = other.content[i][j];
            }
        }
    }
}

// TODO: handle if column does not exist, otherwise it crashes.
// If you pass COL_NAME, you have one of the following situations:
// you ask the original relation, you get exact unique column.
// TODO: add if you ask a derived relation, you may get any match
//  for a column with same name for any of the relations.
//
// If you use full name i.e. R1*R2.COL_NAME
// You ask a reduced relation i.e. R1, you get no match.
// you ask original, you get exact unique column.
// you ask derived, you get any match from inner relations
// TODO: if reduced name exists, should we make sure the parent table is matching?
//  -
Column &Relation::operator[](const std::string &columnName) {
    std::string str_2 = columnName;

    // If it does not exists somewhere,
    // Assume "TABLE_NAME.[COL_NAME]" or "TABLE_NAME.COL_NAME" is passed
    // then try to remove the table name and check for "COL_Name"
    if ((schema.find(str_2) == schema.end())) {

        auto pos = str_2.find('.');

        if (pos != std::string::npos) {
            std::string str_1 = columnName.substr(0, pos);
            if (name == str_1) {
                str_2 = columnName.substr(pos + 1);
                return (*this)[str_2];
            } else {
                str_2 = "...NONE...";
                leaderLogging("Column name does not exist: " + columnName + "\n", true);
            }
        } else {
            str_2 = "...NONE...";
            leaderLogging("Column name does not exist: " + columnName + "\n", true);
        }
    }

    if (isBColumnName(str_2)) {
        return getColumnB(str_2);
    } else {
        return getColumnA(str_2);
    }
}

Column &Relation::operator[](const std::shared_ptr<Expression> columnExpression) {
//    return (*this)[columnExpression->getIdentifier()];
    return (*this)[((AttributeExpr*) &*columnExpression)->getIdentifier()];
}

// BColumn name has no "["
bool Relation::isBColumnName(const std::string &columnName) const {
    // int pos = columnName.find('.');
    // if (pos == -1) {
    //     return columnName[0] == '[' && columnName[columnName.size() - 1] == ']';
    // } else {
    //     return columnName[pos + 1] == '[' && columnName[columnName.size() - 1] == ']';
    // }
    int pos = columnName.find('[');
    return pos != -1;
}

Column &Relation::getColumnA(const std::string &columnName) {
    auto other_ = reinterpret_cast<DerivedColumn<Share, ShareTable> *>(schema[columnName]);
    other_->allocateRowShares(&shareTable,true);
    return *schema[columnName];
}

Column &Relation::getColumnB(const std::string &columnName) {
    auto other_ = reinterpret_cast<DerivedColumn<Share, ShareTable> *>(schema[columnName]);
    other_->allocateRowShares(&shareTable, true);
    return *schema[columnName];
}

std::string Relation::getName() const {
    return name;
}

int Relation::getSize() const {
    return shareTable.numRows;
}

std::vector<std::string> Relation::getACols() const {
    std::vector<std::string> cols = getCols();
    std::vector<std::string> cols_a;

    for (const std::string &col : cols) {
        if (!isBColumnName(col)) {
            cols_a.push_back(col);
        }
    }

    return cols_a;
}

std::vector<std::string> Relation::getBCols() const {
    std::vector<std::string> cols = getCols();
    std::vector<std::string> cols_b;

    for (const std::string &col : cols) {
        if (isBColumnName(col)) {
            cols_b.push_back(col);
        }
    }

    return cols_b;
}

std::vector<std::string> Relation::getCols(bool omitTableName) const {
    std::vector<std::string> cols(schema.size());
    if (omitTableName)
        for (auto it = schema.begin(); it != schema.end(); ++it) 
            cols[it->second->getIndex()] = it->first;
    else 
        for (auto it = schema.begin(); it != schema.end(); ++it) 
            cols[it->second->getIndex()] = name + "." + it->first;
    return cols;
}

/**
 * Initialize column at index with bool secret shares of 1
*/
void Relation::initBColumnOnes(int index) {
    for (int i=0; i<this->getSize(); i++) {
        this->shareTable.content[i][index * CR_P.PARTY_REPLICATION] = 1;
        this->shareTable.content[i][index * CR_P.PARTY_REPLICATION + 1] = 1;
    }
}

void Relation::initAColumnZeros(int index) {
    for (int i=0; i<this->getSize(); i++) {
        this->shareTable.content[i][index * CR_P.PARTY_REPLICATION] = 0;
        this->shareTable.content[i][index * CR_P.PARTY_REPLICATION + 1] = 0;
    }
}

void Relation::initAColumnOnes(int index) {
    for (int i=0; i<this->getSize(); i++) {
        if (get_rank() == 0) {
            this->shareTable.content[i][index * CR_P.PARTY_REPLICATION] = 1;
            this->shareTable.content[i][index * CR_P.PARTY_REPLICATION+1] = 0;
        }
        else if (get_rank() == 1) {
            this->shareTable.content[i][index * CR_P.PARTY_REPLICATION] = 0;
            this->shareTable.content[i][index * CR_P.PARTY_REPLICATION+1] = 0;
        }
        else {
            this->shareTable.content[i][index * CR_P.PARTY_REPLICATION] = 0;
            this->shareTable.content[i][index * CR_P.PARTY_REPLICATION+1] = 1;
        }
    }
}

// TODO: remove the protocol replication from inside this function.
//  - make protocol replication only when passing to basic secrecy core API
std::pair<unsigned int *, int> Relation::getIndices(const std::vector<std::shared_ptr<Expression>> &colExpressions) {
    std::pair<unsigned int *, int> indices_pair;

    // Count cols
    std::vector<Column *> cols;
    for (auto expr: colExpressions) {
        Column &c = (*this)[expr];

        if (c.getIndex() >= 0) {
            cols.push_back(&c);
        }
    }

    // Allocate memory
    indices_pair.second = cols.size();
    indices_pair.first = new unsigned int[indices_pair.second];

    // Assign values
    for (int i = 0; i < indices_pair.second; ++i) {
        indices_pair.first[i] = cols[i]->getIndex() * CR_P.PARTY_REPLICATION;
    }

    return indices_pair;
}

std::shared_ptr<Data *> Relation::openRelation() const {
    auto data = allocateColumnWiseContent<Data>(shareTable.numRows, shareTable.numCols, 1);
    auto data_ptr = data.get();

    for (auto it = schema.begin(); it != schema.end(); ++it) {
        int colIndex = it->second->getIndex();
        int colIndex_2 = 2 * colIndex;

        if (isBColumnName(it->first)) {
            for (int i = 0; i < shareTable.numRows; ++i) {
                data_ptr[colIndex][i] = open_b(shareTable.content[i][colIndex_2]);
            }
        } else {
            for (int i = 0; i < shareTable.numRows; ++i) {
                data_ptr[colIndex][i] = open_a(shareTable.content[i][colIndex_2]);
            }
        }
    }
    return data;
}
