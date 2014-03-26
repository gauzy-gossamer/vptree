/*
 *  vptree.h
 *  Implementation of a vantage-point tree.
 *
 *  Created by Laurens van der Maaten.
 *  Copyright 2012, Delft University of Technology. All rights reserved.
 *
 */


#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <queue>
#include <limits>
#include <float.h>


#ifndef VPTREE_H
#define VPTREE_H

int debug = 0;

class DataPoint
{
    int _D;
    int _ind;
    double* _x;

public:
    DataPoint() {
        _D = 1;
        _ind = -1;
        _x = NULL;
    }
    DataPoint(int D, int ind, double* x) {
        _D = D;
        _ind = ind;
        _x = (double*) malloc(_D * sizeof(double));
        for(int d = 0; d < _D; d++) _x[d] = x[d];
    }
    DataPoint(const DataPoint& other) {                     // this makes a deep copy -- should not free anything
        if(this != &other) {
            _D = other.dimensionality();
            _ind = other.index();
            _x = (double*) malloc(_D * sizeof(double));      
            for(int d = 0; d < _D; d++) _x[d] = other.x(d);
        }
    }
    ~DataPoint() { if(_x != NULL) free(_x); }
    DataPoint& operator= (const DataPoint& other) {         // asignment should free old object
        if(this != &other) {
            if(_x != NULL) free(_x);
            _D = other.dimensionality();
            _ind = other.index();
            _x = (double*) malloc(_D * sizeof(double));
            for(int d = 0; d < _D; d++) _x[d] = other.x(d);
        }
        return *this;
    }
    int index() const { return _ind; }
    int dimensionality() const { return _D; }
    double x(int d) const { return _x[d]; }

    void view() {
        for (int d = 0; d < dimensionality(); d++) {
            printf("v %i\t%.2f\n", d, this->x(d));
        }
    }

    double distance (const DataPoint &t2) const {
        double dd = .0;
        for(int d = 0; d < dimensionality(); d++) {
            dd += (this->x(d) - t2.x(d)) * (this->x(d) - t2.x(d));
        }
        return sqrt(dd);
    }
};

typedef struct {
    int index;
    double value;
} vals_t;

class DataPointSparse
{
    int _D;
    int _ind;
    vals_t *_x;


public:
    DataPointSparse() {
        _D = 1;
        _ind = -1;
        _x = NULL;
    }

    DataPointSparse(int D, int ind, double *x) {
        _D = D;
        _ind = ind;
        _x = (vals_t*) calloc(_D, sizeof(vals_t));
        int j = 0;
        for(int d = 0; d < _D*2; d += 2) {
            if (debug) {
                printf("ins %i\t%.4f\n", (int)x[d], x[d + 1]);
            }
            _x[j].index = (int)x[d];
            _x[j].value = x[d + 1];
            j++;
        }
    }

    ~DataPointSparse() {
        if (debug) {
            printf("release DataPointSparse %i\n", _D);
        }
        if (_x != NULL) { free(_x); }
    }

    DataPointSparse(const DataPointSparse& other) {
        if(this != &other) {
            _D = other.dimensionality();
            _ind = other.index();
            _x = (vals_t*) calloc(_D, sizeof(vals_t));      
            vals_t *_x_other = other.get_vals_ref();

            for(int d = 0; d < _D; d++) {
                _x[d].index = _x_other[d].index;
                _x[d].value = _x_other[d].value;
            }
        }
    }

    DataPointSparse& operator= (const DataPointSparse& other) {
        if (debug) {
            printf("assign sparse\n");
        }
        if(this != &other) {
            if(_x != NULL) {free(_x); }
            _D = other.dimensionality();
            _ind = other.index();
            _x = (vals_t*) calloc(_D, sizeof(vals_t));
            vals_t *_x_other = other.get_vals_ref();

            for(int d = 0; d < _D; d++) {
                _x[d].index = _x_other[d].index;
                _x[d].value = _x_other[d].value;
            }
        }
        return *this;
    }

    vals_t *get_vals_ref() const { return _x; }

    int dimensionality() const { return _D; }
    int index() const { return _ind; }

    double distance(const DataPointSparse &t2) const {
        int t1_size = this->dimensionality() - 1, t2_size = t2.dimensionality() - 1;
        double dd = 0.0;
        int i = 0, j = 0;
        vals_t *vals = this->get_vals_ref(), *vals2 = t2.get_vals_ref();

        if (debug) {
            printf("dims = %i %i\n", t1_size, t2_size);
        }

        while (true) {
            if (vals[i].index == vals2[j].index) {
                dd += (vals[i].value - vals2[j].value) * (vals[i].value - vals2[j].value);

                if (debug) {
                    printf("tt %i d=%.2f\n", vals2[j].index, dd);
                }
                i++;
                j++;
            } else if (vals[i].index > vals2[j].index) {
                dd += (vals2[j].value)*(vals2[j].value);
                if (debug) {
                    printf("t2 %i %.2f d=%.2f\n", vals2[j].index, vals2[j].value, dd);
                }
                j++;
            } else {
                dd += (vals[i].value)*(vals[i].value);

                if (debug) {
                    printf("t1 %i %.2f d=%.2f\n", vals[i].index, vals[i].value, dd);
                }
                i++;
            }

            if (i > t1_size || j > t2_size) {
                break;
            }
        }
        while (i < t1_size) {
            dd += vals[i].value*vals[i].value;
            if (debug) {
                printf("post t1 %i %.2f d=%.2f\n", vals[i].index, vals[i].value,dd);
            }
            i++;
        }
        while (j < t2_size) {
            dd += vals2[j].value*vals2[j].value;
            if (debug) {
                printf("post t2 %i %.2f d=%.2f\n", vals2[j].index, vals2[j].value, dd);
            }
            j++;
        }

        return sqrt(dd);
    }

};

template<typename T> //, double (*distance)( const T&, const T& )>
class VpTree
{
public:
    
    // Default constructor
    VpTree() : _root(0), sparse(false) {}
    
    // Destructor
    ~VpTree() {
        delete _root;
        //delete _items;
    }

    bool is_sparse() {
        return sparse;
    }

    void set_sparse() {
        sparse = true;
    }

    void view(double *vals, int N) {
        int i, j = 0;
        for (i = 0; i < N*2; i += 2) {
            printf("%.2f\t%.2f\n", vals[i], vals[i+1]);
            j++;
        }
    }

    void create(double *items, int D, int N) {
        this->_D = D;
        std::vector<T> *obj_X = new std::vector<T>(N, T(D, -1, items));
        for(int n = 0; n < N; n++) {
            obj_X->assign(n, T(D, n, items + n * D));
        }

        obj_X->at(0).view();


        this->create(obj_X);
    }

    // Function to create a new VpTree from data
    void create(std::vector<T>* items) {
        delete _root;
        //if (_items) { 
        //    delete _items;
       // }
        _items = items;
        this->_D = _items->at(0).dimensionality();
        if (debug) {
            printf("dim %i\n", this->_D);
        }
        _root = buildFromPoints(0, items->size());
    }

    /* indices instead of numbers */
    /* X is a vector of features. k - number of neighbors to look for */
    void search(double *X, int k, int D, double epsilon, std::vector<int> *indices, std::vector<double> *distances) {
        //if (D != this->_D) {
        //    fprintf(stderr, "incorrect number of dimensions, original was %i\n", this->_D);
        //    return;
        //}
        if (_root == 0) {
            return;
        }
        std::vector<T> results;
        if (debug) {
            printf("search %i d %i\n", k, D);
        }
        T item = T(D, -1, X);

        this->search(item, k, epsilon, &results, distances);

        for(typename std::vector<T>::iterator it = results.begin(); it != results.end(); ++it) {
            indices->push_back(it->index());
        }
    }
    
    // Function that uses the tree to find the k nearest neighbors of target
    void search(const T& target, int k, double epsilon, std::vector<T>* results, std::vector<double>* distances)
    {
        if (_root == 0) {
            return;
        }
        // Use a priority queue to store intermediate results on
        std::priority_queue<HeapItem> heap;
        
        // Variable that tracks the distance to the farthest point in our results
        _tau = DBL_MAX;

        // Perform the search
        search(_root, target, k, epsilon, heap);
        
        // Gather final results
        results->clear(); distances->clear();
        while(!heap.empty()) {
            results->push_back(_items->at(heap.top().index));
            distances->push_back(heap.top().dist);
            heap.pop();
        }
        
        // Results are in reverse order
        std::reverse(results->begin(), results->end());
        std::reverse(distances->begin(), distances->end());
    }
    
private:
    std::vector<T> *_items;
    double _tau;
    int _D;
    bool sparse;
    
    // Single node of a VP tree (has a point and radius; left children are closer to point than the radius)
    struct Node
    {
        int index;              // index of point in node
        double threshold;       // radius(?)
        Node* left;             // points closer by than threshold
        Node* right;            // points farther away than threshold
        
        Node() :
        index(0), threshold(0.), left(0), right(0) {}
        
        ~Node() {               // destructor
            delete left;
            delete right;
        }
    }* _root;
    
    
    // An item on the intermediate result queue
    struct HeapItem {
        HeapItem( int index, double dist) :
        index(index), dist(dist) {}
        int index;
        double dist;
        bool operator<(const HeapItem& o) const {
            return dist < o.dist;
        }
    };
    
    // Distance comparator for use in std::nth_element
    struct DistanceComparator
    {
        const T& item;
        DistanceComparator(const T& item) : item(item) {}
        bool operator()(const T& a, const T& b) {
            return item.distance(a) < item.distance(b);
        }
    };
    
    // Function that (recursively) fills the tree
    Node* buildFromPoints( int lower, int upper )
    {
        if (upper == lower) {     // indicates that we're done here!
            return NULL;
        }
        
        // Lower index is center of current node
        Node* node = new Node();
        node->index = lower;
        
        // if we did not arrive at leaf yet
        if (upper - lower > 1) {
            // Choose an arbitrary point and move it to the start
            int i = (int) ((double)rand() / RAND_MAX * (upper - lower - 1)) + lower;
            //printf("rand %i\n", i);

            int median = (upper + lower) / 2;
            //printf("dist %.2f; %i; %i\n", _items->at(i).distance(_items->at(median)), i, median);

            std::swap(_items->at(lower), _items->at(i));

            
            // Partition around the median distance
            std::nth_element(_items->begin() + lower + 1,
                             _items->begin() + median,
                             _items->begin() + upper,
                             DistanceComparator(_items->at(lower)));
            
            // Threshold of the new node will be the distance to the median
            node->threshold = _items->at(lower).distance(_items->at(median));
            
            // Recursively build tree
            node->index = lower;
            node->left = buildFromPoints(lower + 1, median);
            node->right = buildFromPoints(median, upper);
        }
        
        // Return result
        return node;
    }
    
    // Helper function that searches the tree    
    void search(Node* node, const T& target, int k, double epsilon, std::priority_queue<HeapItem>& heap)
    {
        if(node == NULL) return;     // indicates that we're done here
        
        // Compute distance between target and current node
        double dist = target.distance(_items->at(node->index));

        // If current node within radius tau
        if(dist < _tau && dist <= epsilon) {
            // remove furthest node from result list (if we already have k results)
            if(heap.size() == k) heap.pop();
            // add current node to result list
            heap.push(HeapItem(node->index, dist));
            // update value of tau (farthest point in result list)
            if(heap.size() == k) _tau = heap.top().dist;
        }
        
        // Return if we arrived at a leaf
        if(node->left == NULL && node->right == NULL) {
            return;
        }
        
        // If the target lies within the radius of ball
        if(dist < node->threshold) {
            // if there can still be neighbors inside the ball, recursively search left child first
            if(dist - _tau <= node->threshold) {
                search(node->left, target, k, epsilon, heap);
            }

            // if there can still be neighbors outside the ball, recursively search right child
            if(dist + _tau >= node->threshold) {
                search(node->right, target, k, epsilon, heap);
            }
        
        // If the target lies outsize the radius of the ball
        } else {
            // if there can still be neighbors outside the ball, recursively search right child first
            if(dist + _tau >= node->threshold) {
                search(node->right, target, k, epsilon, heap);
            }
            
            // if there can still be neighbors inside the ball, recursively search left child
            if (dist - _tau <= node->threshold) {
                search(node->left, target, k, epsilon, heap);
            }
        }
    }
};
 
#endif
