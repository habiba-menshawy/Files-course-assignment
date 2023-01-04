#include <iostream>
#include <iostream>
#include <fstream>
#include <cmath>
#include "cstring"
#include "algorithm"
#include "bits/stdc++.h"
using namespace std;
class block {
public:
    int nodeNum;
    int *record;
    int recNum;

    block() {}

    block(int m, int nextAvailable, int recNum) {
        this->nodeNum = m * 2;
        this->recNum = recNum;
        record = new int[(m * 2) + 2];
        record[0] = 0;
        record[1] = -1;
        record[2] = nextAvailable;
        for (int i = 3; i < (nodeNum + 2); i++) {
            record[i] = -1;
        }
    }

    void display() {
        for (int i = 1; i < nodeNum + 2; i++) {
            cout << record[i] << ",";

        }
        cout << endl;
    }

};

void CreateIndexFileFile(char *filename, int numberOfRecords, int m) {
    fstream BtreeIndex;
    BtreeIndex.open(filename, ios::in | ios::binary | ios::out | ios::trunc);
    if (!BtreeIndex.is_open()) cout << "Error creating file" << endl;
    else {
        BtreeIndex.seekp(0, ios::beg);
        for (int i = 0; i < numberOfRecords - 1; i++) {
            block rec(m, i + 1, numberOfRecords);
            BtreeIndex.write((char *) &rec, sizeof(block));
        }
        block rec(m, -1, numberOfRecords);
        BtreeIndex.write((char *) &rec, sizeof(block));
        BtreeIndex.close();
    }

}

void DisplayIndexFileContent(char *filename) {
    fstream BtreeIndex;
    BtreeIndex.open(filename, ios::in | ios::out | ios::binary | ios::app);
    if (!BtreeIndex.is_open()) cout << "Error creating file" << endl;
    else {
        BtreeIndex.seekg(0, ios::beg);
        block rec;
        while (BtreeIndex.read((char *) &rec, sizeof(block))) {
            rec.display();
        }
    }
    BtreeIndex.close();
}

int SearchARecord(char *filename, int RecordID) {
    fstream BtreeIndex;
    BtreeIndex.open(filename, ios::in | ios::out | ios::binary | ios::app);
    if (!BtreeIndex.is_open()) cout << "Error creating file" << endl;
    else {

        block root;
        BtreeIndex.seekg(sizeof(block));
        BtreeIndex.read((char *) &root, sizeof(block));

        int target = -1, i;
        bool flag = false;
        bool f = false;
        block rec;
        for (i = 2; i <= root.record[0] * 2; i += 2) {
            if (root.record[i] >= RecordID) {
                flag = true;
                target = root.record[i + 1];
                BtreeIndex.seekg(target * sizeof(root), ios::beg);
                BtreeIndex.read((char *) &rec, sizeof(root));
                while (rec.record[1]) {
                    int j;
                    for (j = 2; i <= rec.record[0] * 2; j += 2) {
                        if (rec.record[j] >= RecordID) {
                            target = rec.record[j + 1];
                            BtreeIndex.seekg(target * sizeof(root), ios::beg);
                            BtreeIndex.read((char *) &rec, sizeof(root));
                            break;
                        }
                    }
                }
                break;
            }


        }
        if (flag) {

            int j;
            int rr = rec.record[0] * 2;
            for (j = 2; j <= rec.record[0] * 2; j += 2) {
                if (rec.record[j] == RecordID) {
                    target = rec.record[j + 1];
                    f = true;
                    break;
                }
            }

        }
        if (f)
            return target;
        else
            return -1;
    }

}

void sortBlock(block &rec) {
    map<int, int> nodes;
    for (int i = 2; i <= rec.record[0] * 2; i += 2) {
        int r = i + 1;
        nodes[rec.record[i]] = rec.record[r];
    }
    int j = 2;
    for (auto i = nodes.begin(); i != nodes.end(); ++i) {
        rec.record[j] = i->first;
        j = j + 1;
        rec.record[j] = i->second;
        j = j + 1;

    }
}

pair<bool, vector<int>> searchInsert(block &root, fstream &BtreeIndex, int RecordID) {
    vector<int> visited;
    visited.push_back(1);
    int target = -1, i;
    bool flag = false;
    for (i = 2; i <= root.record[0] * 2 && root.record[1] == 1; i += 2) {
        if (root.record[i] >= RecordID) {
            flag = true;
            target = root.record[i + 1];
            visited.push_back(target);
            BtreeIndex.seekg(target * sizeof(root), ios::beg);
            block rec;
            BtreeIndex.read((char *) &rec, sizeof(root));
            while (rec.record[1]) {
                int j;
                for (j = 2; j <= rec.record[0] * 2; j += 2) {
                    if (rec.record[j] >= RecordID) {
                        target = rec.record[j + 1];
                        visited.push_back(target);
                        BtreeIndex.seekg(target * sizeof(root), ios::beg);
                        BtreeIndex.read((char *) &rec, sizeof(root));
                        break;
                    }
                }
            }
            break;
        }

    }
    if (!flag && root.record[1] == 1) {
        target = root.record[--i];
        visited.push_back(target);
        block n;
        BtreeIndex.seekg(target * sizeof(root), ios::beg);
        BtreeIndex.read((char *) &n, sizeof(root));
        while (n.record[1]) {
//            int j;
//            for (j = 0; j <= n.record[0] * 2; j += 2) {
//                target = n.record[j + 1];
//            }
            target=n.record[n.record[0]*2+1];
            visited.push_back(target);
            BtreeIndex.seekg(target * sizeof(root), ios::beg);
            BtreeIndex.read((char *) &n, sizeof(root));

        }
    }
    return {flag, visited};

}

int updateNextAvailableInsert(fstream &BtreeIndex) {
    block header;
    BtreeIndex.seekg(0);
    BtreeIndex.read((char *) &header, sizeof(block));
    int track = header.record[2];
    block check;
    BtreeIndex.seekg(track * sizeof(block));
    BtreeIndex.read((char *) &check, sizeof(block));
    int leaf = check.record[2];
    header.record[2] = leaf;
    return leaf;
}

void updateNextAvailableDelete(fstream &BtreeIndex, int n, block &rec) {
    block header;
    BtreeIndex.seekg(0);
    BtreeIndex.read((char *) &header, sizeof(block));
    if (n < header.record[2] || header.record[2] == -1) {
        rec.record[2] = header.record[2];
        header.record[2] = n;
    } else {
        int track = header.record[2];
        while (true) {
            block check;
            BtreeIndex.seekg(track * sizeof(block));
            BtreeIndex.read((char *) &check, sizeof(block));
            int leaf = check.record[1];
            if (leaf == -1) {
                if (n < check.record[2] || check.record[2] == -1) {
                    rec.record[2] = check.record[2];
                    check.record[2] = n;
                    return;
                }
                track = check.record[2];
            }
        }
    }
}

void updateGreatest(vector<int> visited, fstream &BtreeIndex, int RecordID) {
    block rec;
    int track = visited[visited.size() - 2];
    BtreeIndex.seekg(sizeof(block) * track);
    BtreeIndex.read((char *) &rec, sizeof(block));
    int nn = rec.record[0];
    rec.record[nn * 2] = RecordID;
    if (track == 1)
        return;
    int c = 3;
    while (true) {
        track = visited[visited.size() - c];
        BtreeIndex.seekg(sizeof(block) * track);
        BtreeIndex.read((char *) &rec, sizeof(block));
        int nn = rec.record[0];
        rec.record[nn * 2] = RecordID;
        c++;
        if (track == 1)
            break;
    }
}

void updateDeleted(vector<int> visited, fstream &BtreeIndex, int oldId, int newId) {
    block rec;
    int track = visited[visited.size() - 2];
    BtreeIndex.seekg(sizeof(block) * track);
    BtreeIndex.read((char *) &rec, sizeof(block));
    int nn = rec.record[0];
    int k;
    int target;
    bool flag = false;
    for (k = 2; k <= rec.record[0] * 2; k += 2) {
        int f = rec.record[k];
        if (rec.record[k] == oldId) {
            flag = true;
            target = rec.record[k + 1];
            break;
        }
    }
    if (flag) { rec.record[k] = newId; }
    else {
        return;
    }
    flag = false;
    if (track == 1)
        return;
    int c = 3;
    while (true) {
        int track = visited[visited.size() - c];
        BtreeIndex.seekg(sizeof(block) * track);
        BtreeIndex.read((char *) &rec, sizeof(block));
        for (k = 2; k <= rec.record[0] * 2; k += 2) {
            if (rec.record[k] == oldId) {
                flag = true;
                target = rec.record[k + 1];
                break;
            }
        }
        if (flag) { rec.record[k] = newId; }
        else {
            return;
        }
        flag = false;
        c++;
        if (track == 1)
            break;
    }
}

void insertInTo(block &rec, int RecordID, int Reference) {
    rec.record[0] = rec.record[0] + 1;
    int nn = rec.record[0];
    rec.record[nn * 2] = RecordID;
    rec.record[nn * 2 + 1] = Reference;
    sortBlock(rec);
}


int splitRoot(block &rec, fstream &BtreeIndex, int RecordID, int Reference, int leaf, block *temp1 = nullptr,
              block *temp2 = nullptr, bool isNotGreatest = false, int next = 0) {
    block header;
    BtreeIndex.seekg(0);
    BtreeIndex.read((char *) &header, sizeof(block));
    int nextAvailable = header.record[2];
    if (nextAvailable == -1) {
        cout << "Can't be entered" << endl;
        return -1;
    } else {

        block rec1;
        block rec2;
        int next1 = nextAvailable;
        BtreeIndex.seekg(nextAvailable * sizeof(block), ios::beg);
        BtreeIndex.read((char *) &rec1, sizeof(block));
        if (nextAvailable == -1) {
            cout << "Can't be entered" << endl;
            return -1;
        }
        nextAvailable = updateNextAvailableInsert(BtreeIndex);
        if (nextAvailable == -1) {
            cout << "Can't be entered" << endl;
            return -1;
        }
        int next2 = nextAvailable;
        BtreeIndex.seekg(nextAvailable * sizeof(block), ios::beg);
        BtreeIndex.read((char *) &rec2, sizeof(block));
        nextAvailable = updateNextAvailableInsert(BtreeIndex);

        insertInTo(rec, RecordID, Reference);
        rec1.record[1] = leaf;
        rec2.record[1] = leaf;

        if ((rec.nodeNum / 2) % 2 == 0) {
            int nn = (rec.nodeNum / 2) + 4;
            rec1.record[0] = floor(rec.record[0] / 2);
            rec2.record[0] = floor(rec.record[0] / 2)-1;
            for (int i = 2; i < nn; i++) {
                rec1.record[i] = rec.record[i];
            }
            int j = 2;
            for (int i = nn; i < rec.nodeNum + 4; i++) {
                rec2.record[j] = rec.record[i];
                j = j + 1;
            }

            rec.record[0] = 2;
            rec.record[1] = 1;

            rec.record[2] = rec1.record[nn - 2];
            rec.record[3] = next1;

            rec.record[4] = rec2.record[nn - 4];
            rec.record[5] = next2;
        } else {
            int nn = (rec.nodeNum / 2) + 3;
            rec1.record[0] = floor(rec.record[0] / 2);
            rec2.record[0] = floor(rec.record[0] / 2);
            for (int i = 2; i < nn; i++) {
                rec1.record[i] = rec.record[i];
            }
            int j = 2;
            for (int i = nn; i < rec.nodeNum + 4; i++) {
                rec2.record[j] = rec.record[i];
                j = j + 1;
            }

            rec.record[0] = 2;
            rec.record[1] = 1;

            rec.record[2] = rec1.record[nn - 2];
            rec.record[3] = next1;

            rec.record[4] = rec2.record[nn - 2];
            rec.record[5] = next2;
        }
        for (int i = 6; i < rec.nodeNum + 2; i++) {
            rec.record[i] = -1;
        }

        if (leaf == 1) {
            map<int, int> rec1mp;
            for (int i = 2; i <= rec1.record[0] * 2; i += 2) {
                int r = i + 1;
                rec1mp[rec1.record[i]] = rec1.record[r];
            }
            map<int, int> rec2mp;
            for (int i = 2; i <= rec2.record[0] * 2; i += 2) {
                int r = i + 1;
                rec2mp[rec2.record[i]] = rec2.record[r];
            }

            int tmp1n = temp1->record[0];
            int tmp2n = temp2->record[0];
            if (rec1mp.find(temp1->record[tmp1n * 2]) == rec1mp.end() &&
                rec1mp.find(temp2->record[tmp2n * 2]) == rec1mp.end()) {
                // not found

                if (isNotGreatest) {
                    rec2mp[temp1->record[tmp1n * 2]] = rec2mp[temp2->record[tmp2n * 2]];
                    rec2mp[temp2->record[tmp2n * 2]] = next;
                } else {

                    rec2mp[temp1->record[tmp1n * 2]] = rec2mp[temp2->record[tmp2n * 2 - 2]];
                    rec2mp.erase(temp2->record[tmp2n * 2 - 2]);
                    rec2mp[temp2->record[tmp2n * 2]] = next;

                }

                int j = 2;
                for (auto i = rec2mp.begin(); i != rec2mp.end(); ++i) {
                    rec2.record[j] = i->first;
                    j = j + 1;
                    rec2.record[j] = i->second;
                    j = j + 1;

                }
            } else if (rec2mp.find(temp1->record[tmp1n * 2]) == rec2mp.end() &&
                       rec2mp.find(temp2->record[tmp2n * 2]) == rec2mp.end()) {
                // found

                if (isNotGreatest) {
                    rec1mp[temp1->record[tmp1n * 2]] = rec1mp[temp2->record[tmp2n * 2]];
                    rec1mp[temp2->record[tmp2n * 2]] = next;
                } else {
                    rec1mp[temp1->record[tmp1n * 2]] = rec1mp[temp2->record[tmp2n * 2 - 2]];
                    rec1mp.erase(rec1mp[temp2->record[tmp2n * 2 - 2]]);
                    rec1mp[temp2->record[tmp2n * 2]] = next;
                }

                int j = 2;
                for (auto i = rec1mp.begin(); i != rec1mp.end(); ++i) {
                    rec1.record[j] = i->first;
                    j = j + 1;
                    rec1.record[j] = i->second;
                    j = j + 1;

                }
            } else {
                //1 fe 7eta w el tany fe 7eta
                /*/////////////////////////////////////////////////////*/
                //rec1 3aizeen 7agat temp1 bs
                if (isNotGreatest) {
                    rec1mp[temp1->record[tmp1n * 2]] = rec2mp[temp2->record[tmp2n * 2]];
                    rec1mp.erase(temp2->record[tmp2n * 2]);

                } else {
                    rec1mp[temp1->record[tmp1n * 2]] = rec2mp[temp2->record[tmp2n * 2 - 2]];


                }
                if (isNotGreatest) {

                    rec2mp[temp2->record[tmp2n * 2]] = next;
                } else {

                    rec2mp.erase(temp2->record[tmp2n * 2 - 2]);
                    rec2mp[temp2->record[tmp2n * 2]] = next;

                }
                int j = 2;
                for (auto i = rec2mp.begin(); i != rec2mp.end(); ++i) {
                    rec2.record[j] = i->first;
                    j = j + 1;
                    rec2.record[j] = i->second;
                    j = j + 1;

                }
                j = 2;
                for (auto i = rec1mp.begin(); i != rec1mp.end(); ++i) {
                    rec1.record[j] = i->first;
                    j = j + 1;
                    rec1.record[j] = i->second;
                    j = j + 1;

                }

            }
        }
    }
}

pair<block, block> splitNorm(block &rec, fstream &BtreeIndex, int RecordID, int Reference) {
    block header;
    BtreeIndex.seekg(0);
    BtreeIndex.read((char *) &header, sizeof(block));
    int nextAvailable = header.record[2];

    if (nextAvailable == -1) {
        cout << "Can't be entered" << endl;
    } else {
        updateNextAvailableInsert(BtreeIndex);//////////////////////////
        block root;
        BtreeIndex.seekg(sizeof(header));
        BtreeIndex.read((char *) &root, sizeof(block));
        insertInTo(rec, RecordID, Reference);
        block rec1;
        BtreeIndex.seekg(nextAvailable * sizeof(block), ios::beg);
        BtreeIndex.read((char *) &rec1, sizeof(block));
        rec1.record[1] = 0;
        int nn;
        if ((rec.nodeNum / 2) % 2 == 0) {
            nn = (rec.nodeNum / 2) + 4;
            rec1.record[0] = floor(rec.record[0] / 2) - 1;
            rec.record[0] = floor(rec.record[0] / 2);

            int j = 2;
            for (int i = nn; i < rec.nodeNum + 4; i++) {
                rec1.record[j] = rec.record[i];
                j = j + 1;
            }

        } else {
            nn = (rec.nodeNum / 2) + 3;
            rec1.record[0] = floor(rec.record[0] / 2);
            rec.record[0] = floor(rec.record[0] / 2);

            int j = 2;
            for (int i = nn; i < rec.nodeNum + 4; i++) {
                rec1.record[j] = rec.record[i];
                j = j + 1;
            }

        }
        for (int i = nn; i < rec.nodeNum + 2; i++) {
            rec.record[i] = -1;
        }
        return {rec, rec1};

    }
}

pair<block, block> splitParent(block &parent, fstream &BtreeIndex, block temp1, block temp2, bool isNotGreatest,int next) {
    block header;
    BtreeIndex.seekg(0);
    BtreeIndex.read((char *) &header, sizeof(block));
    int nextAvailable = header.record[2];
    if (nextAvailable == -1) {
        cout << "Can't be entered" << endl;
    } else {
        block root;
        BtreeIndex.seekg(sizeof(header));
        BtreeIndex.read((char *) &root, sizeof(block));

        map<int, int> Nmp;
        for (int i = 2; i <= parent.record[0] * 2; i += 2) {
            int r = i + 1;
            Nmp[parent.record[i]] = parent.record[r];
        }
        int rec1n = temp1.record[0];
        int rec2n = temp2.record[0];
        // found
            if (isNotGreatest) {
                Nmp[temp1.record[rec1n * 2]] = Nmp[temp2.record[rec1n * 2]];
                Nmp[temp2.record[rec2n * 2]] = next;
            } else {
                Nmp[temp1.record[rec1n * 2]] = Nmp[temp2.record[rec1n * 2 - 2]];
                Nmp.erase(temp2.record[rec1n * 2 - 2]);
                Nmp[temp2.record[rec2n * 2]] = next;
            }
        int j = 2;
        for (auto i = Nmp.begin(); i != Nmp.end(); ++i) {
            parent.record[j] = i->first;
            j = j + 1;
            parent.record[j] = i->second;
            j = j + 1;

        }
        parent.record[0] = parent.record[0] + 1;
/*//////////////////////////////////*/
        block rec1;
        BtreeIndex.seekg(nextAvailable * sizeof(block), ios::beg);
        BtreeIndex.read((char *) &rec1, sizeof(block));
        nextAvailable = updateNextAvailableInsert(BtreeIndex);
        rec1.record[1] = 0;
        int nn;
        if ((parent.nodeNum / 2) % 2 == 0) {
            nn = (parent.nodeNum / 2) + 4;
            rec1.record[0] = floor(parent.record[0] / 2) - 1;
            parent.record[0] = floor(parent.record[0] / 2);

            int j = 2;
            for (int i = nn; i < parent.nodeNum + 4; i++) {
                rec1.record[j] = parent.record[i];
                j = j + 1;
            }

        } else {
            nn = (parent.nodeNum / 2) + 3;
            rec1.record[0] = floor(parent.record[0] / 2);
            parent.record[0] = floor(parent.record[0] / 2);

            int j = 2;
            for (int i = nn; i < parent.nodeNum + 4; i++) {
                rec1.record[j] = parent.record[i];
                j = j + 1;
            }

        }
        for (int i = nn; i < parent.nodeNum + 2; i++) {
            parent.record[i] = -1;
        }

        return {parent, rec1};
    }
}

int secondaryInsert(fstream &BtreeIndex, int RecordID, int Reference) {
    block header;
    BtreeIndex.seekg(0);
    BtreeIndex.read((char *) &header, sizeof(block));
    int nextAvailable = header.record[2];
    block root;
    BtreeIndex.seekg(sizeof(block));
    BtreeIndex.read((char *) &root, sizeof(block));
    bool isNotGreatest = searchInsert(root, BtreeIndex, RecordID).first;
    vector<int> visited = searchInsert(root, BtreeIndex, RecordID).second;
    int track = visited[visited.size() - 1];
    block rec;
    BtreeIndex.seekg(sizeof(block) * track);
    BtreeIndex.read((char *) &rec, sizeof(block));
    int c = 2;
    if (rec.record[0] < rec.nodeNum / 2) {
        insertInTo(rec, RecordID, Reference);
        if (!isNotGreatest) {
            updateGreatest(visited, BtreeIndex, RecordID);
        }
    } else {
        block n;
        BtreeIndex.seekg(visited[visited.size() - c] * sizeof(root), ios::beg);
        BtreeIndex.read((char *) &n, sizeof(root));
        if (header.record[2] == -1) {
            cout << "Can't be entered" << endl;
            return -1;
        }

        int next = header.record[2];
        pair<block, block> temps = splitNorm(rec, BtreeIndex, RecordID, Reference);
        block temp1 = temps.first;
        block temp2 = temps.second;

        while (n.record[1]) {
            if (n.record[0] < n.nodeNum / 2) {

                map<int, int> Nmp;
                for (int i = 2; i <= n.record[0] * 2; i += 2) {
                    int r = i + 1;
                    Nmp[n.record[i]] = n.record[r];
                }
                int rec1n = temp1.record[0];
                int rec2n = temp2.record[0];
                // found

                if (isNotGreatest) {
                    Nmp[temp1.record[rec1n * 2]] = Nmp[temp2.record[rec1n * 2]];
                    Nmp[temp2.record[rec2n * 2]] = next;
                } else {
                    Nmp[temp1.record[rec1n * 2]] = Nmp[temp2.record[rec1n * 2 - 2]];
                    Nmp.erase(temp2.record[rec1n * 2 - 2]);
                    Nmp[temp2.record[rec2n * 2]] = next;
                }


                n.record[0] = n.record[0] + 1;
                int j = 2;
                for (auto i = Nmp.begin(); i != Nmp.end(); ++i) {
                    n.record[j] = i->first;
                    j = j + 1;
                    n.record[j] = i->second;
                    j = j + 1;

                }

                if (!isNotGreatest)
                    updateGreatest(visited, BtreeIndex, RecordID);
                break;
            } else {

                if (visited[visited.size() - c] == 1) {
                    if (header.record[2] == -1) {
                        cout << "Can't be entered" << endl;
                        return -1;
                    }
                    splitRoot(root, BtreeIndex, temp1.record[temp1.record[0] * 2], 0, 1, &temp1, &temp2, isNotGreatest,
                              next);
                    if (!isNotGreatest)
                        updateGreatest(visited, BtreeIndex, RecordID);
                    return 1;
                } else {
                    if (header.record[2] == -1) {
                        cout << "Can't be entered" << endl;
                        return -1;
                    }

                    temps = splitParent(n, BtreeIndex, temp1, temp2, isNotGreatest,next);
                    temp1 = temps.first;
                    temp2 = temps.second;
                    BtreeIndex.seekg(visited[visited.size() - c] * sizeof(root), ios::beg);
                    BtreeIndex.read((char *) &n, sizeof(root));
                    next = header.record[2];
                    c++;
                }
            }
        }

    }
}

int InsertNewRecordAtIndex(char *filename, int RecordID, int Reference) {
    fstream BtreeIndex;
    BtreeIndex.open(filename, ios::binary | ios::out | ios::in | ios::app);
    block header;
    BtreeIndex.seekg(0);
    BtreeIndex.read((char *) &header, sizeof(block));
    int nextAvailable = header.record[2];
    block root;
    BtreeIndex.seekg(sizeof(header));
    BtreeIndex.read((char *) &root, sizeof(block));
    if (nextAvailable == 1) {
        //tree in empty
        //NextAvailable
        header.record[2] = header.record[2] + 1;
        //FirstNode
        root.record[2] = RecordID;
        root.record[3] = Reference;
        //LeafSlot
        root.record[1] = 0;
        //counter
        root.record[0] = root.record[0] + 1;
        BtreeIndex.close();
        return 1;
    } else {
        if (root.record[1] == 0) {

            if (root.record[0] < root.nodeNum / 2) {

                insertInTo(root, RecordID, Reference);
                BtreeIndex.close();
                return 1;
            } else {
                splitRoot(root, BtreeIndex, RecordID, Reference, 0);
                vector<int> visited = searchInsert(root, BtreeIndex, RecordID).second;
                BtreeIndex.close();
                return visited[visited.size() - 1];
            }
        } else {
            if (secondaryInsert(BtreeIndex, RecordID, Reference) == -1) {
                return -1;
            }
            vector<int> visited = searchInsert(root, BtreeIndex, RecordID).second;
            BtreeIndex.close();
            return visited[visited.size() - 1];
        }
    }
    BtreeIndex.close();
}

void deleteNode(fstream &BtreeIndex, int RecordID, vector<int> visited, block &rec) {
    int target;
    int j;
    int rr = rec.record[0] * 2;
    for (j = 2; j <= rec.record[0] * 2; j += 2) {
        if (rec.record[j] == RecordID) {
            target = rec.record[j + 1];
            break;
        }
    }
    map<int, int> nodes;
    for (int i = 2; i <= rec.record[0] * 2; i += 2) {
        int r = i + 1;
        nodes[rec.record[i]] = rec.record[r];

    }
    nodes.erase(RecordID);
    if (target == rec.record[rec.record[0] * 2 + 1]) {
        updateDeleted(visited, BtreeIndex, RecordID, rec.record[(rec.record[0] * 2) - 2]);
    }
    rec.record[0] = rec.record[0] - 1;
    j = 2;
    for (auto i = nodes.begin(); i != nodes.end(); ++i) {
        rec.record[j] = i->first;
        j = j + 1;
        rec.record[j] = i->second;
        j = j + 1;

    }
    for (int i = j; i < rec.nodeNum + 2; i++) {
        rec.record[i] = -1;
    }

}

void mergeParents(fstream &BtreeIndex, block &root, vector<int> visited) {
    int c;
    block parent;
    int parentTrack = visited[visited.size() - 2];
    BtreeIndex.seekg(sizeof(block) * parentTrack);
    BtreeIndex.read((char *) &parent, sizeof(block));
    if (parent.record[0] >= parent.nodeNum / 4 || visited[visited.size() - 2] == 1) {
        return;
    }
    block grandParent;
    int grandParentTrack = visited[visited.size() - 3];
    BtreeIndex.seekg(sizeof(block) * grandParentTrack);
    BtreeIndex.read((char *) &grandParent, sizeof(block));
    int j;
    int target;
    bool flag = false;
    for (j = 2; j <= grandParent.record[0] * 2; j += 2) {
        if (grandParent.record[j] == parent.record[parent.record[0] * 2]) {
            flag = true;
            break;
        }
    }


    if (flag) {
        if (j != 2) {
            pair<bool, vector<int>> sib1Res = searchInsert(root, BtreeIndex, grandParent.record[j - 2]);
            vector<int> visSib1 = sib1Res.second;
            block sib1;
            int sib1Track = visSib1[visSib1.size() - 2];
            BtreeIndex.seekg(sizeof(block) * sib1Track);
            BtreeIndex.read((char *) &sib1, sizeof(block));
            if (sib1.record[0] <= sib1.nodeNum / 4) {
                int y = grandParent.record[grandParent.record[0] * 2 - 2];
                int b = sib1.record[sib1.record[0] * 2];
                int z = parent.record[parent.record[0] * 2];
                int rec1h = 2;
                int temp = visited[visited.size() - 1];


                visited.pop_back();
                for (auto x: visited) {
                    block del;
                    BtreeIndex.seekg(sizeof(block) * x);
                    BtreeIndex.read((char *) &del, sizeof(block));
                    if (x != visited[visited.size() - 1])
                        deleteNode(BtreeIndex, parent.record[parent.record[0] * 2], visited, del);
                }
                visSib1.pop_back();
                updateDeleted(visSib1, BtreeIndex, sib1.record[sib1.record[0] * 2],
                              parent.record[parent.record[0] * 2]);
                visited.emplace_back(temp);
                while (rec1h <= parent.record[0] * 2) {
                    insertInTo(sib1, parent.record[rec1h], parent.record[rec1h + 1]);
                    rec1h = rec1h + 2;
                }
                for (int i = 0; i < parent.recNum + 2; i++) {
                    parent.record[i] = -1;
                }
                parent.record[0] = 0;
                updateNextAvailableDelete(BtreeIndex, parentTrack, parent);
                return;
            }
            if (sib1.record[0] > sib1.nodeNum / 4) {
                //a3ml insert el awl abl ma ams7
                insertInTo(parent, sib1.record[sib1.record[0] * 2], sib1.record[sib1.record[0] * 2 + 1]);
                deleteNode(BtreeIndex, sib1.record[sib1.record[0] * 2], visSib1, sib1);

                return;
            }

        }

        if (j != grandParent.record[0] * 2) {
            pair<bool, vector<int>> sib2Res = searchInsert(root, BtreeIndex, grandParent.record[j + 2]);
            vector<int> visSib2 = sib2Res.second;
            block sib2;
            int sib2Track = grandParent.record[j + 3];
            BtreeIndex.seekg(sizeof(block) * sib2Track);
            BtreeIndex.read((char *) &sib2, sizeof(block));
            if (sib2.record[0] <= sib2.nodeNum / 4) {

                int rec1h = 2;
                while (rec1h <= parent.record[0] * 2) {
                    insertInTo(sib2, parent.record[rec1h], parent.record[rec1h + 1]);
                    rec1h = rec1h + 2;
                }
                for (auto x: visited) {
                    block del;
                    BtreeIndex.seekg(sizeof(block) * x);
                    BtreeIndex.read((char *) &del, sizeof(block));
                    if (x != visited[visited.size() - 1])
                        deleteNode(BtreeIndex, parent.record[parent.record[0] * 2], visited, del);
                }
                for (int i = 0; i < parent.recNum + 2; i++) {
                    parent.record[i] = -1;
                }
                parent.record[0] = 0;
                updateNextAvailableDelete(BtreeIndex, parentTrack, parent);
                return;
            }

            if (sib2.record[0] > sib2.nodeNum / 4) {
                //a3ml insert el awl abl ma ams7
                int temp = visited[visited.size() - 1];
                visited.pop_back();
                insertInTo(parent, sib2.record[2], sib2.record[3]);
                updateDeleted(visited, BtreeIndex, parent.record[parent.record[0] * 2 - 2], sib2.record[2]);
                deleteNode(BtreeIndex, sib2.record[2], visSib2, sib2);
                visited.emplace_back(temp);
                return;
            }

        }
    }
}

void mergeSibling(fstream &BtreeIndex, block &root, int RecordID, vector<int> visited, bool isNotGreatest, block &rec) {
    int trackRec = visited[visited.size() - 1];
    block parent;
    int parentTrack = visited[visited.size() - 2];
    BtreeIndex.seekg(sizeof(block) * parentTrack);
    BtreeIndex.read((char *) &parent, sizeof(block));
    int j;
    int target;
    bool flag = false;
    for (j = 2; j <= parent.record[0] * 2; j += 2) {
        if (parent.record[j] == rec.record[rec.record[0] * 2]) {
            flag = true;
            break;
        }
    }
    if (flag) {
        if (j != 2) {
            pair<bool, vector<int>> sib1Res = searchInsert(root, BtreeIndex, parent.record[j - 2]);
            vector<int> visSib1 = sib1Res.second;
            block sib1;
            int sib1Track = visSib1[visSib1.size() - 1];
            BtreeIndex.seekg(sizeof(block) * sib1Track);
            BtreeIndex.read((char *) &sib1, sizeof(block));
            if (sib1.record[0] <= sib1.nodeNum / 4) {
                for (auto x: visited) {
                    if (x == 1)
                        continue;
                    block del;
                    BtreeIndex.seekg(sizeof(block) * x);
                    BtreeIndex.read((char *) &del, sizeof(block));
                    if (x == visited[visited.size() - 1]) {
                        deleteNode(BtreeIndex, RecordID, visited, del);
                    } else {
                        deleteNode(BtreeIndex, rec.record[rec.record[0] * 2], visited, del);
                    }
                }
                int rec1h = 2;
                updateDeleted(visSib1, BtreeIndex, sib1.record[sib1.record[0] * 2], rec.record[rec.record[0] * 2]);
                updateDeleted(visited, BtreeIndex, sib1.record[sib1.record[0] * 2], rec.record[rec.record[0] * 2]);
                int ex = 0;
                for (int h = sib1.record[0] * 2 + 2; rec1h <= rec.record[0] * 2 + 1; h += 1) {
                    ex++;
                    sib1.record[h] = rec.record[rec1h];
                    rec1h = rec1h + 1;
                }
                sortBlock(sib1);
                sib1.record[0] = sib1.record[0] + ex - 1;

                for (int i = 0; i < rec.recNum + 2; i++) {
                    rec.record[i] = -1;
                }
                rec.record[0] = 0;
                updateNextAvailableDelete(BtreeIndex, trackRec, rec);
                while (visited[visited.size() - 1] != 1) {
                    mergeParents(BtreeIndex, root, visited);
                    visited.pop_back();
                }
                return;
            }
            if (sib1.record[0] > sib1.nodeNum / 4) {
                deleteNode(BtreeIndex, RecordID, visited, rec);
                //a3ml insert el awl abl ma ams7
                insertInTo(rec, sib1.record[sib1.record[0] * 2], sib1.record[sib1.record[0] * 2 + 1]);
                //   updateGreatest(visited, BtreeIndex, sib1.record[sib1.record[0] * 2]);
                deleteNode(BtreeIndex, sib1.record[sib1.record[0] * 2], visSib1, sib1);
                return;
            }
        }

        if (j != parent.record[0] * 2) {
            pair<bool, vector<int>> sib2Res = searchInsert(root, BtreeIndex, parent.record[j + 2]);
            vector<int> visSib2 = sib2Res.second;
            block sib2;
            int sib2Track = visSib2[visSib2.size() - 1];
            BtreeIndex.seekg(sizeof(block) * sib2Track);
            BtreeIndex.read((char *) &sib2, sizeof(block));
            if (sib2.record[0] <= sib2.nodeNum / 4) {
                for (auto x: visited) {
                    if (x == 1)
                        continue;
                    block del;
                    BtreeIndex.seekg(sizeof(block) * x);
                    BtreeIndex.read((char *) &del, sizeof(block));
                    if (x == visited[visited.size() - 1]) {
                        deleteNode(BtreeIndex, RecordID, visited, del);
                    } else {
                        deleteNode(BtreeIndex, rec.record[rec.record[0] * 2], visited, del);
                    }
                }

                int rec1h = 2;
                updateDeleted(visSib2, BtreeIndex, rec.record[rec.record[0] * 2], sib2.record[sib2.record[0] * 2]);
                updateDeleted(visited, BtreeIndex, rec.record[rec.record[0] * 2], sib2.record[sib2.record[0] * 2]);

                while (rec1h <= rec.record[0] * 2) {

                    insertInTo(sib2, rec.record[rec1h], rec.record[rec1h + 1]);
                    rec1h = rec1h + 2;
                }


                for (int i = 0; i < rec.recNum + 2; i++) {
                    rec.record[i] = -1;
                }
                rec.record[0] = 0;
                updateNextAvailableDelete(BtreeIndex, trackRec, rec);

                while (visited[visited.size() - 1] != 1) {
                    mergeParents(BtreeIndex, root, visited);
                    visited.pop_back();
                }
                return;

            }

            if (sib2.record[0] > sib2.nodeNum / 4) {

                deleteNode(BtreeIndex, RecordID, visited, rec);
                //   vector<int> visSib = searchInsert(root, BtreeIndex, sib1.record[sib1.record[0] * 2]).second;
                //a3ml insert el awl abl ma ams7
                insertInTo(rec, sib2.record[2], sib2.record[3]);
                updateDeleted(visited, BtreeIndex, rec.record[rec.record[0] * 2 - 2], sib2.record[2]);
                deleteNode(BtreeIndex, sib2.record[2], visSib2, sib2);
                return;
            }
        }
    }
}

void DeleteRecordFromIndex(char *filename, int RecordID) {
    fstream BtreeIndex;
    BtreeIndex.open(filename, ios::binary | ios::out | ios::in | ios::app);
    block header;
    BtreeIndex.seekg(0);
    BtreeIndex.read((char *) &header, sizeof(block));
    int nextAvailable = header.record[2];
    block root;
    BtreeIndex.seekg(sizeof(header));
    BtreeIndex.read((char *) &root, sizeof(block));

    if (nextAvailable == 1) {
        cout << "Tree is empty No ID to delete" << endl;
        return;
    }
    if (SearchARecord(filename, RecordID) == -1) {
        cout << "This ID doesn't exist in the index" << endl;
        return;
    }
    pair<bool, vector<int>> searchRes = searchInsert(root, BtreeIndex, RecordID);
    bool isNotGreatest = searchRes.first;
    vector<int> visited = searchRes.second;
    block rec;
    int track = visited[visited.size() - 1];
    BtreeIndex.seekg(sizeof(block) * track);
    BtreeIndex.read((char *) &rec, sizeof(block));
    if (track == 1) {
        deleteNode(BtreeIndex, RecordID, visited, rec);
        return;
    }
    if (rec.record[0] == 1) {
        deleteNode(BtreeIndex, RecordID, visited, rec);
        return;
    }
    if (rec.record[0] > rec.nodeNum / 4) {
        deleteNode(BtreeIndex, RecordID, visited, rec);
        return;
    } else {
        block parent;
        int parentTrack = visited[visited.size() - 2];
        BtreeIndex.seekg(sizeof(block) * parentTrack);
        BtreeIndex.read((char *) &parent, sizeof(block));
        mergeSibling(BtreeIndex, root, RecordID, visited, isNotGreatest, rec);

    }

    BtreeIndex.close();
}
int main() {
    CreateIndexFileFile("BtreeIndex.bin", 10, 5);
    // DisplayIndexFileContent("BtreeIndex.bin");
    cout<<endl;
    InsertNewRecordAtIndex("BtreeIndex.bin", 3, 12);
    InsertNewRecordAtIndex("BtreeIndex.bin", 7, 24);
    InsertNewRecordAtIndex("BtreeIndex.bin", 10, 48);
    InsertNewRecordAtIndex("BtreeIndex.bin", 24, 60);
    InsertNewRecordAtIndex("BtreeIndex.bin", 14, 72);
    DisplayIndexFileContent("BtreeIndex.bin");
    cout << endl;
    InsertNewRecordAtIndex("BtreeIndex.bin", 19, 84);
    DisplayIndexFileContent("BtreeIndex.bin");
    cout << endl;
    InsertNewRecordAtIndex("BtreeIndex.bin", 30, 196);
    InsertNewRecordAtIndex("BtreeIndex.bin", 15, 108);
    InsertNewRecordAtIndex("BtreeIndex.bin", 1, 120);
    InsertNewRecordAtIndex("BtreeIndex.bin", 5, 132);
    DisplayIndexFileContent("BtreeIndex.bin");
    cout << endl;
    InsertNewRecordAtIndex("BtreeIndex.bin", 2, 144);
    DisplayIndexFileContent("BtreeIndex.bin");
    cout << endl;
    InsertNewRecordAtIndex("BtreeIndex.bin", 8, 156);
    InsertNewRecordAtIndex("BtreeIndex.bin", 9, 168);
    InsertNewRecordAtIndex("BtreeIndex.bin", 6, 180);
    InsertNewRecordAtIndex("BtreeIndex.bin", 11, 192);
    InsertNewRecordAtIndex("BtreeIndex.bin", 12, 204);
    InsertNewRecordAtIndex("BtreeIndex.bin", 17, 216);
    InsertNewRecordAtIndex("BtreeIndex.bin", 18, 228);
    DisplayIndexFileContent("BtreeIndex.bin");
    cout << endl;
    InsertNewRecordAtIndex("BtreeIndex.bin", 32, 240);
    DisplayIndexFileContent("BtreeIndex.bin");
    cout << endl;
    cout << SearchARecord("BtreeIndex.bin", 3) << endl;
    cout << SearchARecord("BtreeIndex.bin", 7) << endl;
    cout << SearchARecord("BtreeIndex.bin", 10) << endl;
    cout << SearchARecord("BtreeIndex.bin", 24) << endl;
    cout << SearchARecord("BtreeIndex.bin", 14) << endl;
    cout << SearchARecord("BtreeIndex.bin", 19) << endl;
    cout << SearchARecord("BtreeIndex.bin", 30) << endl;
    cout << SearchARecord("BtreeIndex.bin", 15) << endl;
    cout << SearchARecord("BtreeIndex.bin", 1) << endl;
    cout << SearchARecord("BtreeIndex.bin", 5) << endl;
    cout << SearchARecord("BtreeIndex.bin", 2) << endl;
    cout << SearchARecord("BtreeIndex.bin", 8) << endl;
    cout << SearchARecord("BtreeIndex.bin", 9) << endl;
    cout << SearchARecord("BtreeIndex.bin", 6) << endl;
    cout << SearchARecord("BtreeIndex.bin", 11) << endl;
    cout << SearchARecord("BtreeIndex.bin", 12) << endl;
    cout << SearchARecord("BtreeIndex.bin", 17) << endl;
    cout << SearchARecord("BtreeIndex.bin", 18) << endl;
    cout << SearchARecord("BtreeIndex.bin", 32) << endl;
    cout << SearchARecord("BtreeIndex.bin", 100) << endl;
    cout << SearchARecord("BtreeIndex.bin", 0) << endl;
    cout << SearchARecord("BtreeIndex.bin", 31) << endl;
    cout << endl;
    DeleteRecordFromIndex("BtreeIndex.bin", 10);
    DisplayIndexFileContent("BtreeIndex.bin");
    cout << endl;
    DeleteRecordFromIndex("BtreeIndex.bin", 9);
    DisplayIndexFileContent("BtreeIndex.bin");
    cout << endl;
    DeleteRecordFromIndex("BtreeIndex.bin", 8);
    DisplayIndexFileContent("BtreeIndex.bin");
    cout << endl;
    return 0;
}
