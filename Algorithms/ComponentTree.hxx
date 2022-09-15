/*
 * This file is part of libTIM.
 *
 * Copyright (©) 2005-2013  Benoit Naegel
 * Copyright (©) 2013 Theo de Carpentier
 * Copyright (©) 2022-2023  Cyril Meyer
 */

#include <functional>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <utility>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace LibTIM {

using std::map;
using std::vector;

template <class T>
ComponentTree<T>::ComponentTree(Image<T>& img) : m_root(0), m_img(img) {
  FlatSE connexity;
  connexity.make2DN8();
  SalembierRecursiveImplementation<T> strategy(this, connexity);

  m_root = strategy.computeTree();
  strategy.computeAttributes(m_root);
}

template <class T>
ComponentTree<T>::ComponentTree(Image<T>& img, FlatSE& connexity)
    : m_root(0), m_img(img) {
  SalembierRecursiveImplementation<T> strategy(this, connexity);

  m_root = strategy.computeTree();
  strategy.computeAttributes(m_root);
}

template <class T>
ComponentTree<T>::ComponentTree(Image<T>& img, FlatSE& connexity,
                                unsigned int delta)
    : m_root(0), m_img(img) {
  SalembierRecursiveImplementation<T> strategy(this, connexity);

  m_root = strategy.computeTree();
  strategy.computeAttributes(m_root, delta);
}

template <class T>
ComponentTree<T>::ComponentTree(Image<T>& img, FlatSE& connexity,
                                ComputedAttributes ca, unsigned int delta)
    : m_root(0),
      m_img(img)

{
  SalembierRecursiveImplementation<T> strategy(this, connexity);

  m_root = strategy.computeTree();

  if (ca & ComputedAttributes::OTSU) {
    computeNeighborhoodAttributes(delta);
  }

  strategy.computeAttributes(m_root, ca, delta);
}

template <class T>
ComponentTree<T>::~ComponentTree() {
  erase_tree();
}

template <class T>
int ComponentTree<T>::computeNeighborhoodAttributes(int r) {
  FlatSE se;
  se.make2DEuclidianBall(r);

  std::queue<Node*> fifo;
  fifo.push(m_root);

  Image<bool> active(m_img.getSize());
  active.fill(0);

  // for each node, without recursive call
  while (!fifo.empty()) {
    // n = current node
    Node* n = fifo.front();
    fifo.pop();

    // active pixels are neighborhs pixels
    active.fill(true);
    std::vector<TOffset> pixels = merge_pixels(n);
    for (int i = 0; i < pixels.size(); i++) active(pixels[i]) = false;

    vector<int> ndg;

    for (int i = 0; i < pixels.size(); i++) {
      Point<TCoord> p = m_img.getCoord(pixels[i]);

      for (int j = 0; j < se.getNbPoints(); j++) {
        Point<TCoord> q = p + se.getPoint(j);

        if (m_img.isPosValid(q)) {
          if (active(q) == true) {
            n->area_nghb += 1;
            n->sum_nghb += m_img(q);
            n->sum_square_nghb += m_img(q) * m_img(q);

            ndg.push_back(m_img(q));

            active(q) = false;
          }
        }
      }
    }

    if (n->area_nghb > 0) {
      n->mean_nghb = (long double)n->sum_nghb / (long double)n->area_nghb;

      n->variance_nghb =
          ((long double)n->sum_square_nghb / (long double)n->area_nghb) -
          n->mean_nghb * n->mean_nghb;
    }

    for (int i = 0; i < n->childs.size(); i++) {
      fifo.push(n->childs[i]);
    }
  }

  return 0;
}

template <class T>
void ComponentTree<T>::erase_tree() {
  int tot = 0;
  if (m_root != 0) {
    std::queue<Node*> fifo;
    fifo.push(m_root);
    while (!fifo.empty()) {
      Node* curNode = fifo.front();
      fifo.pop();

      for (int i = 0; i < curNode->childs.size(); i++) {
        fifo.push(curNode->childs[i]);
      }
      delete curNode;
      tot++;
    }
  }
}

template <class T>
Image<T>& ComponentTree<T>::constructImageOptimized() {
  int numberNonActives = 0;
  if (m_root != 0)
    if (m_root->active == true) {
      std::queue<Node*> fifo;
      std::queue<Node*> fifoChilds;
      fifo.push(m_root);
      while (!fifo.empty()) {
        Node* tmp = fifo.front();
        fifo.pop();

        for (std::vector<Node*>::iterator it = tmp->childs.begin();
             it != tmp->childs.end(); ++it) {
          if ((*it)->active == false) {
            fifoChilds.push(*it);

            while (!fifoChilds.empty()) {
              Node* child = fifoChilds.front();
              fifoChilds.pop();

              for (std::vector<TOffset>::iterator it3 = child->pixels.begin();
                   it3 != child->pixels.end(); ++it3)
                m_img(*it3) = (T)tmp->h;

              for (std::vector<Node*>::iterator it2 = child->childs.begin();
                   it2 != child->childs.end(); ++it2) {
                fifoChilds.push(*it2);
              }
            }

          } else
            fifo.push(*it);
        }
      }
    } else
      m_img.fill(T(0));
  else
    m_img.fill(T(0));

  return m_img;
}

template <class T>
void ComponentTree<T>::constructImageMin(Image<T>& res) {
  if (m_root->active == true) {
    std::queue<Node*> fifo;
    fifo.push(m_root);
    while (!fifo.empty()) {
      Node* tmp = fifo.front();
      fifo.pop();

      // update pixels contained in the node
      for (std::vector<TOffset>::iterator it = tmp->pixels.begin();
           it != tmp->pixels.end(); ++it)
        res(*it) = (T)tmp->h;

      for (std::vector<Node*>::iterator it = tmp->childs.begin();
           it != tmp->childs.end(); ++it) {
        // if child->active is false, "cut" the subtree and hence search
        // all pixels of all subnodes
        if ((*it)->active == false) {
          std::vector<TOffset> subPixels = merge_pixels(*it);
          for (std::vector<TOffset>::iterator it = subPixels.begin();
               it != subPixels.end(); ++it)
            res(*it) = (T)tmp->h;
        } else
          fifo.push(*it);
      }
    }
  } else
    res.fill(T(0));
}

// Does not work!!!!!
// TODO: implement construct image Max
template <class T>
void ComponentTree<T>::constructImageMax(Image<T>& res) {
  res.fill(0);

  std::queue<Node*> fifoLeafs;
  std::queue<Node*> fifo;
  fifo.push(m_root);

  while (!fifo.empty()) {
    Node* currentNode = fifo.front();
    fifo.pop();

    currentNode->status = true;
    if (currentNode->childs.size() != 0) {
      for (Node::ContainerChilds::iterator it = currentNode->childs.begin();
           it != currentNode->childs.end(); ++it)
        fifo.push(*it);
    } else
      fifoLeafs.push(currentNode);
  }
  while (!fifoLeafs.empty()) {
    Node* tmp = fifoLeafs.front();
    fifoLeafs.pop();

    if (tmp->active == false && tmp->father->status == true) {
      // we put father in the queue
      fifoLeafs.push(tmp->father);
      tmp->father->status = false;
    } else if (tmp->active == true) {
      // we construct image node at this level
      std::vector<TOffset> subPixels = merge_pixels(tmp);
      for (std::vector<TOffset>::iterator it = subPixels.begin();
           it != subPixels.end(); ++it)
        res(*it) = (T)tmp->h;
    }
  }
}

template <class T>
void ComponentTree<T>::constructImageDirect(Image<T>& res) {
  res.fill(T(0));

  std::queue<Node*> fifo;
  fifo.push(m_root);
  while (!fifo.empty()) {
    Node* tmp = fifo.front();
    fifo.pop();

    if (tmp->active == true) {
      for (std::vector<TOffset>::iterator it = tmp->pixels.begin();
           it != tmp->pixels.end(); ++it)
        res(*it) = (T)tmp->h;

      for (std::vector<Node*>::iterator it = tmp->childs.begin();
           it != tmp->childs.end(); ++it) {
        // return all pixels of all consecutive false subnodes
        // stop when an active node is found
        if ((*it)->active == false) {
          std::vector<TOffset> subPixels = merge_pixelsFalseNodes(*it);
          for (std::vector<TOffset>::iterator it2 = subPixels.begin();
               it2 != subPixels.end(); ++it2)
            res(*it2) = (T)tmp->h;
        }

        fifo.push(*it);
      }
    } else {
      for (std::vector<Node*>::iterator it = tmp->childs.begin();
           it != tmp->childs.end(); ++it) {
        fifo.push(*it);
      }
    }
  }
}

template <class T>
void ComponentTree<T>::constructImageDirectExpe(Image<T>& res) {
  res.fill(0);

  std::queue<Node*> fifo;
  std::queue<Node*> fifo2;

  fifo.push(m_root);

  while (!fifo.empty()) {
    Node* tmp = fifo.front();
    fifo.pop();

    if (tmp->active == true) {
      fifo2.push(tmp);
    }

    else
      for (std::vector<Node*>::iterator it = tmp->childs.begin();
           it != tmp->childs.end(); ++it) {
        fifo.push(*it);
      }
  }

  while (!fifo2.empty()) {
    Node* tmp = fifo2.front();
    fifo2.pop();

    for (std::vector<TOffset>::iterator it = tmp->pixels.begin();
         it != tmp->pixels.end(); ++it)
      res(*it) = (T)tmp->h;

    for (std::vector<Node*>::iterator it = tmp->childs.begin();
         it != tmp->childs.end(); ++it) {
      if ((*it)->active == false) {
        (*it)->h = tmp->h;
      }

      fifo2.push(*it);
    }
  }
}

template <class T>
Image<T> ComponentTree<T>::constructImage(ConstructionDecision decision) {
  Image<T> res(m_img.getSize());

  if (m_root != 0) {
    switch (decision) {
      case MIN:
        constructImageMin(res);
        break;
      case MAX:
        constructImageMax(res);
        break;
      case DIRECT:
        constructImageDirectExpe(res);
        break;
    }
  }

  else
    res.fill(T(0));

  return res;
}

template <class T>
template <class TVal>
TVal ComponentTree<T>::getAttribute(Node* n,
                                    ComponentTree::Attribute attribute_id) {
  switch (attribute_id) {
    case H:
      return n->h;
    case AREA:
      return n->area;
    case AREA_D_AREAN_H:
      return n->area_derivative_areaN_h;
    case AREA_D_AREAN_H_D:
      return n->area_derivative_areaN_h_derivative;
    case AREA_D_H:
      return n->area_derivative_h;
    case AREA_D_AREAN:
      return n->area_derivative_areaN;
    case MSER:
      return n->mser;
    case AREA_D_DELTA_H:
      return n->area_derivative_delta_h;
    case AREA_D_DELTA_AREAF:
      return n->area_derivative_delta_areaF;
    case MEAN:
      return n->mean;
    case VARIANCE:
      return n->variance;
    case MEAN_NGHB:
      return n->mean_nghb;
    case VARIANCE_NGHB:
      return n->variance_nghb;
    case OTSU:
      return n->otsu;
    case CONTRAST:
      return n->contrast;
    case VOLUME:
      return n->volume;
    case MGB:
      return n->mean_gradient_border;
    case CONTOUR_LENGTH:
      return n->contourLength;
    case COMPLEXITY:
      return n->complexity;
    case COMPACITY:
      return n->compacity;
  }
  return 0;
}

template <class T>
template <class TVal, class TSel>
void ComponentTree<T>::constructImageAttributeMin(
    Image<TVal>& res, ComponentTree::Attribute value_attribute,
    ComponentTree::Attribute selection_attribute) {
  res.fill(TVal(0));

  std::vector<Node*> nodes = indexedNodes();
  for (TSize i = 0; i < res.getSizeX(); i++)
    for (TSize j = 0; j < res.getSizeY(); j++)
      for (TSize k = 0; k < res.getSizeZ(); k++) {
        Node* n = indexedCoordToNode(i, j, k, nodes);
        // noeud selectionné
        Node* n_s = n;
        TSel attr = getAttribute<TSel>(n, selection_attribute);
        // minimum, dans la branche parent
        TSel attr_father;
        // parcours de l'arbre
        while (n->father != m_root) {
          n = n->father;
          attr_father = getAttribute<TSel>(n, selection_attribute);

          if (attr_father < attr && attr_father > 0) {
            n_s = n;
            attr = attr_father;
          }
        }

        res(i, j, k) = getAttribute<TVal>(n_s, value_attribute);
      }
}

template <class T>
template <class TVal, class TSel>
void ComponentTree<T>::constructImageAttributeMax(
    Image<TVal>& res, ComponentTree::Attribute value_attribute,
    ComponentTree::Attribute selection_attribute) {
  res.fill(TVal(0));

  std::vector<Node*> nodes = indexedNodes();
  for (TSize i = 0; i < res.getSizeX(); i++)
    for (TSize j = 0; j < res.getSizeY(); j++)
      for (TSize k = 0; k < res.getSizeZ(); k++) {
        Node* n = indexedCoordToNode(i, j, k, nodes);
        // noeud selectionné
        Node* n_s = n;
        TSel attr = getAttribute<TSel>(n, selection_attribute);
        // maximum dans la branche parent
        TSel attr_father;
        // parcours de l'arbre
        while (n->father != m_root) {
          n = n->father;
          attr_father = getAttribute<TSel>(n, selection_attribute);

          if (attr_father > attr &&
              attr_father < std::numeric_limits<TSel>::max()) {
            n_s = n;
            attr = attr_father;
          }
        }

        res(i, j, k) = getAttribute<TVal>(n_s, value_attribute);
      }
}

template <class T>
template <class TVal>
void ComponentTree<T>::constructImageAttributeDirect(
    Image<TVal>& res, ComponentTree::Attribute value_attribute) {
  res.fill(TVal(0));

  std::vector<Node*> nodes = indexedNodes();
  for (TSize i = 0; i < res.getSizeX(); i++)
    for (TSize j = 0; j < res.getSizeY(); j++)
      for (TSize k = 0; k < res.getSizeZ(); k++) {
        Node* n = indexedCoordToNode(i, j, k, nodes);
        TVal attr = getAttribute<TVal>(n, value_attribute);
        res(i, j, k) = attr;
      }
}

template <class T>
template <class TVal, class TSel, class TLimit>
void ComponentTree<T>::constructImageAttributeMin(
    Image<TVal>& res, ComponentTree::Attribute value_attribute,
    ComponentTree::Attribute selection_attribute, Attribute limit_attribute,
    TLimit limit_min, TLimit limit_max) {
  res.fill(TVal(0));

  std::vector<Node*> nodes = indexedNodes();
  for (TSize i = 0; i < res.getSizeX(); i++)
    for (TSize j = 0; j < res.getSizeY(); j++)
      for (TSize k = 0; k < res.getSizeZ(); k++) {
        Node* n = indexedCoordToNode(i, j, k, nodes);
        // limit min
        while (n->father != m_root &&
               getAttribute<TLimit>(n->father, limit_attribute) < limit_min) {
          n = n->father;
        }
        // noeud selectionné
        Node* n_s = n;
        TSel attr = getAttribute<TSel>(n, selection_attribute);
        // minimum, dans la branche parent
        TSel attr_father;
        // parcours de l'arbre et limit max
        while (n->father != m_root &&
               getAttribute<TLimit>(n->father, limit_attribute) < limit_max) {
          n = n->father;
          attr_father = getAttribute<TSel>(n, selection_attribute);

          if (attr_father < attr && attr_father > 0) {
            n_s = n;
            attr = attr_father;
          }
        }

        res(i, j, k) = getAttribute<TVal>(n_s, value_attribute);
      }
}

template <class T>
template <class TVal, class TSel, class TLimit>
void ComponentTree<T>::constructImageAttributeMax(
    Image<TVal>& res, ComponentTree::Attribute value_attribute,
    ComponentTree::Attribute selection_attribute, Attribute limit_attribute,
    TLimit limit_min, TLimit limit_max) {
  res.fill(TVal(0));

  std::vector<Node*> nodes = indexedNodes();
  for (TSize i = 0; i < res.getSizeX(); i++)
    for (TSize j = 0; j < res.getSizeY(); j++)
      for (TSize k = 0; k < res.getSizeZ(); k++) {
        Node* n = indexedCoordToNode(i, j, k, nodes);
        // limit min
        while (n->father != m_root &&
               getAttribute<TLimit>(n->father, limit_attribute) < limit_min) {
          n = n->father;
        }
        // noeud selectionné
        Node* n_s = n;
        TSel attr = getAttribute<TSel>(n, selection_attribute);
        // maximum dans la branche parent
        TSel attr_father;
        // parcours de l'arbre et limit max
        while (n->father != m_root &&
               getAttribute<TLimit>(n->father, limit_attribute) < limit_max) {
          n = n->father;
          attr_father = getAttribute<TSel>(n, selection_attribute);

          if (attr_father > attr &&
              attr_father < std::numeric_limits<TSel>::max()) {
            n_s = n;
            attr = attr_father;
          }
        }

        res(i, j, k) = getAttribute<TVal>(n_s, value_attribute);
      }
}

template <class T>
template <class TVal, class TLimit>
void ComponentTree<T>::constructImageAttributeDirect(
    Image<TVal>& res, ComponentTree::Attribute value_attribute,
    Attribute limit_attribute, TLimit limit_min, TLimit limit_max) {
  res.fill(TVal(0));

  std::vector<Node*> nodes = indexedNodes();
  for (TSize i = 0; i < res.getSizeX(); i++)
    for (TSize j = 0; j < res.getSizeY(); j++)
      for (TSize k = 0; k < res.getSizeZ(); k++) {
        Node* n = indexedCoordToNode(i, j, k, nodes);
        // limit min
        while (n->father != m_root &&
               getAttribute<TLimit>(n->father, limit_attribute) < limit_min) {
          n = n->father;
        }
        TVal attr = getAttribute<TVal>(n, value_attribute);
        res(i, j, k) = attr;
      }
}

template <class T>
template <class TVal, class TSel>
Image<TVal> ComponentTree<T>::constructImageAttribute(
    ComponentTree::Attribute value_attribute,
    ComponentTree::Attribute selection_attribute,
    ComponentTree::ConstructionDecision selection_rule) {
  Image<TVal> res(m_img.getSize());

  if (m_root != 0) {
    switch (selection_rule) {
      case MIN:
        constructImageAttributeMin<TVal, TSel>(res, value_attribute,
                                               selection_attribute);
        break;
      case MAX:
        constructImageAttributeMax<TVal, TSel>(res, value_attribute,
                                               selection_attribute);
        break;
      case DIRECT:
        constructImageAttributeDirect<TVal>(res, value_attribute);
        break;
    }
  } else
    res.fill(TVal(0));

  return res;
}

template <class T>
template <class TVal, class TSel, class TLimit>
Image<TVal> ComponentTree<T>::constructImageAttribute(
    ComponentTree::Attribute value_attribute,
    ComponentTree::Attribute selection_attribute,
    ComponentTree::ConstructionDecision selection_rule,
    Attribute limit_attribute, TLimit limit_min, TLimit limit_max) {
  Image<TVal> res(m_img.getSize());

  if (m_root != 0) {
    switch (selection_rule) {
      case MIN:
        constructImageAttributeMin<TVal, TSel, TLimit>(
            res, value_attribute, selection_attribute, limit_attribute,
            limit_min, limit_max);
        break;
      case MAX:
        constructImageAttributeMax<TVal, TSel, TLimit>(
            res, value_attribute, selection_attribute, limit_attribute,
            limit_min, limit_max);
        break;
      case DIRECT:
        constructImageAttributeDirect<TVal, TLimit>(
            res, value_attribute, limit_attribute, limit_min, limit_max);
        break;
    }
  } else
    res.fill(TVal(0));

  return res;
}

template <class T>
void ComponentTree<T>::constructNode(Image<T>& res, Node* node) {
  std::queue<Node*> fifo;
  fifo.push(node);

  while (!fifo.empty()) {
    Node* tmp = fifo.front();
    fifo.pop();

    for (std::vector<TOffset>::iterator it = tmp->pixels.begin();
         it != tmp->pixels.end(); ++it)
      res(*it) = (T)tmp->h;
    for (std::vector<Node*>::iterator it = tmp->childs.begin();
         it != tmp->childs.end(); ++it)
      fifo.push(*it);
  }
}

template <class T>
void ComponentTree<T>::constructNodeDirect(Image<T>& res, Node* node) {
  std::queue<Node*> fifo;
  fifo.push(node);

  int h = node->h;

  while (!fifo.empty()) {
    Node* tmp = fifo.front();
    fifo.pop();

    for (std::vector<TOffset>::iterator it = tmp->pixels.begin();
         it != tmp->pixels.end(); ++it)
      res(*it) = (T)h;
    for (std::vector<Node*>::iterator it = tmp->childs.begin();
         it != tmp->childs.end(); ++it)
      fifo.push(*it);
  }
}

template <class T>
void ComponentTree<T>::setFalse() {
  if (m_root != 0) {
    std::queue<Node*> fifo;
    fifo.push(m_root);

    while (!fifo.empty()) {
      Node* tmp = fifo.front();
      fifo.pop();

      tmp->active = false;

      for (std::vector<Node*>::iterator it = tmp->childs.begin();
           it != tmp->childs.end(); ++it)
        fifo.push(*it);
    }
  }
}

// Test whether the se is include in the component (pixels)

template <class T>
bool ComponentTree<T>::isInclude(FlatSE& se, Node::ContainerPixels& pixels) {
  // Case where the se is larger than the component:
  // obviously se does not fit in
  if (se.getNbPoints() > pixels.size()) {
    return false;
  } else {
    FlatSE::iterator itSe;
    Node::ContainerPixels::iterator itPixels;
    Node::ContainerPixels::iterator itPixels2;

    for (itPixels = pixels.begin(); itPixels != pixels.end(); ++itPixels) {
      bool isInclude = true;
      for (itSe = se.begin(); itSe != se.end(); ++itSe) {
        TOffset q = *itSe + *itPixels;
        // Test if q belongs to pixels set
        for (itPixels2 = pixels.begin(); itPixels2 != pixels.end(); ++itPixels2)
          if (*itPixels2 == q) break;
        if (itPixels2 == pixels.end()) {
          isInclude = false;
          break;
        }
      }
      if (isInclude == true) return true;
    }

    return false;
  }
}

// aggregate and return all pixels belonging to subtree
template <class T>
std::vector<TOffset> ComponentTree<T>::merge_pixels(Node* tree) {
  vector<TOffset> res;

  std::queue<Node*> fifo;
  fifo.push(tree);
  while (!fifo.empty()) {
    Node* tmp = fifo.front();
    fifo.pop();
    for (int i = 0; i < tmp->pixels.size(); i++) res.push_back(tmp->pixels[i]);
    for (int i = 0; i < tmp->childs.size(); i++) {
      fifo.push(tmp->childs[i]);
    }
  }

  return res;
}

// return all the pixels of the subtree
// stop when it reaches an active node
template <class T>
std::vector<TOffset> ComponentTree<T>::merge_pixelsFalseNodes(Node* tree) {
  vector<TOffset> res;

  std::queue<Node*> fifo;
  fifo.push(tree);
  while (!fifo.empty()) {
    Node* tmp = fifo.front();
    fifo.pop();
    if (tmp->active == false) {
      for (int i = 0; i < tmp->pixels.size(); i++)
        res.push_back(tmp->pixels[i]);
      for (int i = 0; i < tmp->childs.size(); i++) {
        fifo.push(tmp->childs[i]);
      }
    }
  }

  return res;
}

// aggregate and return all pixels belonging to subtree
template <class T>
void ComponentTree<T>::merge_pixels(Node* tree, std::vector<TOffset>& res) {
  std::queue<Node*> fifo;
  fifo.push(tree);
  while (!fifo.empty()) {
    Node* tmp = fifo.front();
    fifo.pop();
    for (int i = 0; i < tmp->pixels.size(); i++) res.push_back(tmp->pixels[i]);
    for (int i = 0; i < tmp->childs.size(); i++) {
      fifo.push(tmp->childs[i]);
    }
  }
}

//////////////////////////////////////////////////////////////
//
//
//
// Component Tree Filters
//
//
//
//////////////////////////////////////////////////////////////

// reinitialization of tree

template <class T>
int ComponentTree<T>::restore() {
  if (m_root != 0) {
    std::queue<Node*> fifo;
    fifo.push(m_root);

    while (!fifo.empty()) {
      Node* curNode = fifo.front();
      fifo.pop();
      curNode->active = true;
      curNode->h = curNode->ori_h;
      Node::ContainerChilds::iterator it;
      for (it = curNode->childs.begin(); it != curNode->childs.end(); ++it)
        fifo.push(*it);
    }
  }
}

template <class T>
int ComponentTree<T>::areaFiltering(int64_t tMin, int64_t tMax) {
  if (m_root != 0) {
    std::queue<Node*> fifo;
    fifo.push(m_root);
    while (!fifo.empty()) {
      Node* curNode = fifo.front();
      fifo.pop();

      if (curNode->area < tMin || curNode->area > tMax) {
        curNode->active = false;
      }

      for (int i = 0; i < curNode->childs.size(); i++) {
        fifo.push(curNode->childs[i]);
      }
    }
  }
  return 0;
}

template <class T>
int ComponentTree<T>::volumicFiltering(int tMin, int tMax) {
  if (m_root != 0) {
    std::queue<Node*> fifo;
    fifo.push(m_root);
    while (!fifo.empty()) {
      Node* curNode = fifo.front();
      fifo.pop();

      if (curNode->volume < tMin || curNode->volume > tMax) {
        curNode->active = false;
      }

      for (int i = 0; i < curNode->childs.size(); i++) {
        fifo.push(curNode->childs[i]);
      }
    }
  }
}

template <class T>
int ComponentTree<T>::contrastFiltering(int tMin, int tMax) {
  if (m_root != 0) {
    std::queue<Node*> fifo;
    fifo.push(m_root);

    while (!fifo.empty()) {
      Node* curNode = fifo.front();
      fifo.pop();
      if (curNode->contrast < tMin || curNode->contrast > tMax) {
        curNode->active = false;
      }

      for (int i = 0; i < curNode->childs.size(); i++) {
        fifo.push(curNode->childs[i]);
      }
    }
  }
  return 0;
}

template <class T>
Node* ComponentTree<T>::coordToNode(TCoord x, TCoord y) {
  TOffset offset = m_img.getOffset(x, y);
  return offsetToNode(offset);
}

template <class T>
Node* ComponentTree<T>::coordToNode(TCoord x, TCoord y, TCoord z) {
  TOffset offset = m_img.getOffset(x, y, z);
  return offsetToNode(offset);
}

template <class T>
Node* ComponentTree<T>::indexedCoordToNode(TCoord x, TCoord y, TCoord z,
                                           std::vector<Node*>& nodes) {
  TOffset offset = m_img.getOffset(x, y, z);
  return nodes[offset];
}

template <class T>
std::vector<Node*> ComponentTree<T>::indexedNodes() {
  unsigned int img_size =
      m_img.getSizeX() * m_img.getSizeY() * m_img.getSizeZ();
  std::vector<Node*> index(img_size);

  std::queue<Node*> fifo;
  fifo.push(m_root);

  while (!fifo.empty()) {
    Node* n = fifo.front();
    fifo.pop();

    if (n != 0) {
      Node::ContainerPixels::iterator it;
      for (it = n->pixels.begin(); it != n->pixels.end(); ++it) {
        index[*it] = n;
      }
      Node::ContainerChilds::iterator jt;
      for (jt = n->childs.begin(); jt != n->childs.end(); ++jt) fifo.push(*jt);
    }
  }
  return index;
}

template <class T>
Node* ComponentTree<T>::offsetToNode(TOffset offset) {
  Node* res = 0;
  std::queue<Node*> fifo;
  fifo.push(m_root);

  while (!fifo.empty()) {
    Node* n = fifo.front();
    fifo.pop();

    if (n != 0) {
      Node::ContainerPixels::iterator it;
      for (it = n->pixels.begin(); it != n->pixels.end(); ++it) {
        if (*it == offset) return n;
      }
      Node::ContainerChilds::iterator jt;
      for (jt = n->childs.begin(); jt != n->childs.end(); ++jt) fifo.push(*jt);
    }
  }
  return res;
}

//////////////////////////////////////////////////////////////
//
//
// Attributes computation
//
//
//
//////////////////////////////////////////////////////////////

template <class T>
int SalembierRecursiveImplementation<T>::computeContrast(Node* tree) {
  if (tree != 0) {
    int current_level = tree->h;
    int current_max = 0;
    int current_contrast = 0;
    Node::ContainerChilds::iterator it;
    for (it = tree->childs.begin(); it != tree->childs.end(); ++it) {
      current_contrast = ((*it)->h - current_level) + computeContrast(*it);
      if (current_contrast > current_max) current_max = current_contrast;
    }
    tree->contrast = current_max;
    return tree->contrast;
  }
  // error
  else
    return -1;
}

template <class T>
int64_t SalembierRecursiveImplementation<T>::computeArea(Node* tree) {
  if (tree != 0) {
    Node::ContainerChilds::iterator it;
    for (it = tree->childs.begin(); it != tree->childs.end(); ++it) {
      tree->area += computeArea(*it);
    }
    return tree->area;
  }
  // error
  else
    return -1;
}

template <class T>
void SalembierRecursiveImplementation<T>::computeAreaDerivative(Node* tree) {
  if (tree != 0) {
    Node::ContainerChilds::iterator it;
    for (it = tree->childs.begin(); it != tree->childs.end(); ++it) {
      computeAreaDerivative(*it);
    }
    tree->area_derivative_areaN_h =
        (((long double)(tree->father->area - tree->area)) /
         ((long double)(tree->h - tree->father->h))) /
        ((long double)(tree->area));
    tree->area_derivative_h = ((long double)(tree->father->area - tree->area)) /
                              ((long double)(tree->h - tree->father->h));
    tree->area_derivative_areaN =
        ((long double)(tree->father->area - tree->area)) /
        ((long double)(tree->area));
  }
}

template <class T>
void SalembierRecursiveImplementation<T>::computeAreaDerivative2(Node* tree) {
  if (tree != 0) {
    Node::ContainerChilds::iterator it;
    for (it = tree->childs.begin(); it != tree->childs.end(); ++it) {
      computeAreaDerivative2(*it);
    }
    tree->area_derivative_areaN_h_derivative =
        tree->father->area_derivative_areaN_h - tree->area_derivative_areaN_h;
  }
}

template <class T>
void SalembierRecursiveImplementation<T>::computeMSER(Node* tree,
                                                      unsigned int delta) {
  if (tree != 0) {
    Node::ContainerChilds::iterator it;
    for (it = tree->childs.begin(); it != tree->childs.end(); ++it) {
      computeMSER(*it, delta);
    }

    tree->mser = std::numeric_limits<long double>::max();
    tree->area_derivative_delta_h = std::numeric_limits<long double>::max();
    tree->area_derivative_delta_areaF = std::numeric_limits<long double>::max();

    Node node = *tree;

    int64_t area_node, area_father;
    int h_node, h_father;

    area_node = node.area;
    h_node = node.h;

    area_father = node.father->area;
    h_father = node.father->h;

    while ((h_node - node.h < (int)delta) &&
           (node.father != node.father->father)) {
      node = *node.father;
    }

    if ((h_node - node.h) >= (int)delta) {
      area_father = node.area;
      h_father = node.h;

      tree->mser =
          ((long double)(area_father - area_node)) / ((long double)(area_node));
      tree->area_derivative_delta_h = ((long double)(area_father - area_node)) /
                                      ((long double)(h_node - h_father));
      tree->area_derivative_delta_areaF =
          ((long double)(area_father - area_node)) /
          ((long double)(area_father));
    }
  }
}

template <class T>
int64_t SalembierRecursiveImplementation<T>::computeSum(Node* tree) {
  if (tree != 0) {
    Node::ContainerChilds::iterator it;
    for (it = tree->childs.begin(); it != tree->childs.end(); ++it) {
      tree->sum += computeSum(*it);
    }
    return tree->sum;
  }
  // error
  else
    return -1;
}

template <class T>
int64_t SalembierRecursiveImplementation<T>::computeSumSquare(Node* tree) {
  if (tree != 0) {
    Node::ContainerChilds::iterator it;
    for (it = tree->childs.begin(); it != tree->childs.end(); ++it) {
      tree->sum_square += computeSumSquare(*it);
    }
    return tree->sum_square;
  }
  // error
  else
    return -1;
}

template <class T>
void SalembierRecursiveImplementation<T>::computeMean(Node* tree) {
  if (tree != 0) {
    Node::ContainerChilds::iterator it;
    for (it = tree->childs.begin(); it != tree->childs.end(); ++it) {
      computeMean(*it);
    }
    tree->mean = (long double)tree->sum / (long double)tree->area;
  }
}

template <class T>
void SalembierRecursiveImplementation<T>::computeVariance(Node* tree) {
  if (tree != 0) {
    Node::ContainerChilds::iterator it;
    for (it = tree->childs.begin(); it != tree->childs.end(); ++it) {
      computeVariance(*it);
    }
    tree->variance = ((long double)tree->sum_square / (long double)tree->area) -
                     tree->mean * tree->mean;
  }
}

template <class T>
void SalembierRecursiveImplementation<T>::computeOtsu(Node* tree) {
  if (tree != 0) {
    Node::ContainerChilds::iterator it;
    for (it = tree->childs.begin(); it != tree->childs.end(); ++it) {
      computeOtsu(*it);
    }
    tree->otsu =
        ((tree->mean - tree->mean_nghb) * (tree->mean - tree->mean_nghb)) /
        (tree->variance + tree->variance_nghb);
  }
}

template <class T>
int64_t SalembierRecursiveImplementation<T>::computeSubNodes(Node* tree) {
  if (tree != 0) {
    Node::ContainerChilds::iterator it;
    for (it = tree->childs.begin(); it != tree->childs.end(); ++it) {
      tree->subNodes = tree->childs.size() + computeSubNodes(*it);
    }
    return tree->subNodes;
  }
  // error
  else
    return -1;
}

/** @brief Compute Volume for each node
 * There is no clear definition of what is volume on discrete image.
 * I choose the following: V(node)=node.area*(node.h-node.father.h)+ sum
 *(V(n.fils)) A special case for the root: the volume is equal to the sum of all
 *grey levels This way the notion is more "intuitive"
 **/

// Warning: this function depends on the attribute "area" of each node
// So the "area" attribute must be computed before any call to this function!
template <class T>
int SalembierRecursiveImplementation<T>::computeVolume(Node* tree) {
  if (tree != 0) {
    int local_contrast = 0;
    // special case for root node
    if (tree->father == tree) {
      local_contrast = tree->h;
    } else
      local_contrast = tree->h - tree->father->h;

    tree->volume = (int)tree->area * local_contrast;

    Node::ContainerChilds::iterator it;
    for (it = tree->childs.begin(); it != tree->childs.end(); ++it) {
      tree->volume += computeVolume(*it);
    }
    return (tree->volume);
  }
  // error
  else
    return -1;
}

template <class T>
void SalembierRecursiveImplementation<T>::computeBorderGradient(Node* tree) {
  if (tree != 0) {
    Node::ContainerChilds::iterator it;
    for (it = tree->childs.begin(); it != tree->childs.end(); ++it) {
      computeBorderGradient(*it);
    }

    long double sum = 0;

    Node::ContainerPixels::iterator itpix;
    for (itpix = tree->pixels_border.begin();
         itpix != tree->pixels_border.end(); ++itpix) {
      sum += imGradient(*itpix);
    }

    tree->mean_gradient_border = sum / tree->pixels_border.size();
  }
}

/** @brief Compute contour length
 *
 **/

// Principle:  -scan image pixels
//             -for each pixel p of value v scan the neighbors (8-connected
//             neighbor gives a 4-connected contour)
//               and retain the neighbor of minimal value v'
//             -update each node n containing p of value v'<vn<v (increment
//             contour length)
//       (if p is not a pixel contour in n, it is not a pixel contour in the
//       ancestors of n)

template <class T>
int SalembierRecursiveImplementation<T>::computeContour(bool save_pixels) {
  // we compute the contour length with STATUS image and m_img
  typename Image<T>::iterator it;
  typename Image<T>::iterator end = imBorder.end();

  this->se.setContext(imBorder.getSize());

  FlatSE::iterator itSE;
  FlatSE::iterator seEnd = se.end();
  ;

  TOffset offset = 0;
  for (it = imBorder.begin(); it != end; ++it, offset++) {
    bool contour = false;
    bool hitsBorder = false;
    T minValue = std::numeric_limits<T>::max();

    if (STATUS(offset) != BORDER_STATUS) {
      for (itSE = se.begin(); itSE != seEnd; ++itSE) {
        TOffset q = offset + *itSE;
        if (STATUS(q) != BORDER_STATUS) {
          if (*it > imBorder(q)) {
            contour = true;
            if (imBorder(q) < minValue) minValue = imBorder(q);
          }
        } else if (STATUS(q) == BORDER_STATUS) {
          // a neighbor of border is a contour point
          // we must propagate this until level Min
          contour = true;
          hitsBorder = true;
          minValue = this->hMin;
        }
      }
    }
    if (contour == true)
    // for all the nodes of levels *it->minValue+1, update contour length
    {
      Node* tmp = index[hToIndex(*it)][STATUS(offset)];
      assert(tmp != 0);
      if (hitsBorder == false)
        while (tmp->h > minValue) {
          tmp->contourLength++;
          if (save_pixels) {
            // conversion offset imBorder->im
            Point<TCoord> imCoord = imBorder.getCoord(offset);
            imCoord.x -= back[0];
            imCoord.y -= back[1];
            imCoord.z -= back[2];
            TOffset imOffset = imCoord.x + imCoord.y * oriSize[0] +
                               imCoord.z * oriSize[0] * oriSize[1];

            tmp->pixels_border.push_back(imOffset);
          }
          tmp = tmp->father;
        }
      else {
        bool stop = false;
        while (!stop) {
          tmp->contourLength++;
          if (save_pixels) {
            // conversion offset imBorder->im
            Point<TCoord> imCoord = imBorder.getCoord(offset);
            imCoord.x -= back[0];
            imCoord.y -= back[1];
            imCoord.z -= back[2];
            TOffset imOffset = imCoord.x + imCoord.y * oriSize[0] +
                               imCoord.z * oriSize[0] * oriSize[1];

            tmp->pixels_border.push_back(imOffset);
          }
          if (tmp != tmp->father)
            tmp = tmp->father;
          else
            stop = true;
        }
      }
    }
  }
  return 1;
}

template <class T>
int SalembierRecursiveImplementation<T>::computeComplexityAndCompacity(
    Node* tree) {
  if (tree != 0) {
    std::queue<Node*> fifo;
    fifo.push(tree);

    while (!fifo.empty()) {
      Node* n = fifo.front();
      fifo.pop();

      if (n->area != 0)
        n->complexity = (int)(1000.0 * n->contourLength / n->area);
      if (n->contourLength != 0) {
        n->compacity = (int)(((double)(4 * M_PI * n->area) /
                              ((double)n->contourLength * n->contourLength)) *
                             1000);
      } else
        n->compacity = 0;
      std::vector<Node*>::iterator it;
      for (it = n->childs.begin(); it != n->childs.end(); ++it) fifo.push(*it);
    }
    return 1;
  }
  // error
  else
    return -1;
}

template <class T>
int SalembierRecursiveImplementation<T>::computeBoundingBox(Node* tree) {
  std::queue<Node*> fifo;
  std::stack<Node*> stackNodes;

  fifo.push(tree);

  while (!fifo.empty()) {
    Node* tmp = fifo.front();
    fifo.pop();

    stackNodes.push(tmp);

    std::vector<Node*>::iterator it;
    for (it = tmp->childs.begin(); it != tmp->childs.end(); ++it)
      fifo.push(*it);
  }

  while (!stackNodes.empty()) {
    Node* tmp = stackNodes.top();
    stackNodes.pop();

    if (tmp->father != tmp) {
      tmp->father->xmin = std::min(tmp->father->xmin, tmp->xmin);
      tmp->father->xmax = std::max(tmp->father->xmax, tmp->xmax);
      tmp->father->ymin = std::min(tmp->father->ymin, tmp->ymin);
      tmp->father->ymax = std::max(tmp->father->ymax, tmp->ymax);
      tmp->father->zmin = std::min(tmp->father->zmin, tmp->zmin);
      tmp->father->zmax = std::max(tmp->father->zmax, tmp->zmax);
    }
  }
  return 1;
}

template <class T>
void SalembierRecursiveImplementation<T>::computeAttributes(Node* tree) {
  if (tree != 0) {
    tree->area = computeArea(tree);
    tree->contrast = computeContrast(tree);
    tree->volume = computeVolume(tree);

    computeContour();
    computeComplexityAndCompacity(tree);
    computeBoundingBox(tree);
    tree->subNodes = computeSubNodes(tree);
  }
}

template <class T>
void SalembierRecursiveImplementation<T>::computeAttributes(
    Node* tree, unsigned int delta) {
  if (tree != 0) {
    tree->area = computeArea(tree);
    computeAreaDerivative(tree);
    computeAreaDerivative2(tree);
    computeMSER(tree, delta);

    tree->contrast = computeContrast(tree);
    tree->volume = computeVolume(tree);
  }
}

template <class T>
void SalembierRecursiveImplementation<T>::computeAttributes(
    Node* tree, ComputedAttributes ca, unsigned int delta) {
  if (tree != 0) {
    if (ca & ComputedAttributes::AREA) {
      tree->area = computeArea(tree);

      if (ca & ComputedAttributes::OTSU) {
        tree->sum = computeSum(tree);
        tree->sum_square = computeSumSquare(tree);
        computeMean(tree);
        computeVariance(tree);
        computeOtsu(tree);
      }
    }
    if (ca & ComputedAttributes::AREA_DERIVATIVES) {
      computeAreaDerivative(tree);
      computeAreaDerivative2(tree);
      computeMSER(tree, delta);
    }
    if (ca & ComputedAttributes::CONTRAST) {
      tree->contrast = computeContrast(tree);
    }
    if (ca & ComputedAttributes::VOLUME) {
      tree->volume = computeVolume(tree);
    }
    if (ca & ComputedAttributes::BORDER_GRADIENT) {
      computeContour(true);
      computeBorderGradient(tree);
    }
    if (ca & ComputedAttributes::COMP_LEXITY_ACITY) {
      if (ca & ComputedAttributes::BORDER_GRADIENT) {
        computeComplexityAndCompacity(tree);
      } else {
        computeContour();
        computeComplexityAndCompacity(tree);
      }
    }
    if (ca & ComputedAttributes::BOUNDING_BOX) {
      computeBoundingBox(tree);
    }
    if (ca & ComputedAttributes::SUB_NODES) {
      tree->subNodes = computeSubNodes(tree);
    }
  }
}

//////////////////////////////////////////////////////////////

template <class T>
inline void SalembierRecursiveImplementation<T>::update_attributes(
    Node* n, TOffset& imBorderOffset) {
  // conversion offset imBorder->im
  Point<TCoord> imCoord = imBorder.getCoord(imBorderOffset);
  imCoord.x -= back[0];
  imCoord.y -= back[1];
  imCoord.z -= back[2];

  TOffset imOffset =
      imCoord.x + imCoord.y * oriSize[0] + imCoord.z * oriSize[0] * oriSize[1];

  n->pixels.push_back(imOffset);
  n->area++;
  n->sum += n->h;
  n->sum_square += (n->h * n->h);

  if (imCoord.x < n->xmin) n->xmin = imCoord.x;
  if (imCoord.x > n->xmax) n->xmax = imCoord.x;

  if (imCoord.y < n->ymin) n->ymin = imCoord.y;
  if (imCoord.y > n->ymax) n->ymax = imCoord.y;

  if (imCoord.z < n->zmin) n->zmin = imCoord.z;
  if (imCoord.z > n->zmax) n->zmax = imCoord.z;
}

template <class T>
inline int SalembierRecursiveImplementation<T>::flood(int h) {
  int m;

  while (!hq[h].empty()) {
    TOffset p = hq[h].front();
    hq[h].pop();

    STATUS(p) = number_nodes[h];

    if (index[h][STATUS(p)] == 0) {
      index[h][STATUS(p)] = this->new_node(indexToH(h), STATUS(p));
      ;
    }

    update_attributes(index[h][STATUS(p)], p);

    FlatSE::iterator it;
    FlatSE::iterator end = se.end();

    for (it = se.begin(); it != end; ++it) {
      TOffset q = p + *it;

      if (STATUS(q) == ACTIVE) {
        hq[hToIndex(imBorder(q))].push(q);
        STATUS(q) = NOT_ACTIVE;

        node_at_level[hToIndex(imBorder(q))] = true;

        if (imBorder(q) > imBorder(p)) {
          m = hToIndex(imBorder(q));

          do {
            m = this->flood(m);

          } while (m != h);
        }
      }
    }
  }
  // End of recursion: we have reached a regional maximum
  number_nodes[h] = number_nodes[h] + 1;

  m = h - 1;
  while (m >= hToIndex(hMin) && node_at_level[m] == false) m--;

  if (m >= hToIndex(hMin)) {
    int i = number_nodes[h] - 1;
    int j = number_nodes[m];
    if (index[m][j] == 0) {
      index[m][j] = new_node(indexToH(m), j);
    }

    this->link_node(index[m][j], index[h][i]);
  } else {
    // The father of root is itself
    index[hToIndex(hMin)][0]->father = index[hToIndex(hMin)][0];
  }
  node_at_level[h] = false;
  return m;
}

template <class T>
Node* SalembierRecursiveImplementation<T>::computeTree() {
  // Put the first pixel with value hMin in the queue
  typename Image<T>::iterator it;
  typename Image<T>::iterator end = imBorder.end();
  TOffset offset = 0;
  for (it = imBorder.begin(); it != end; ++it, offset++)
    if (*it == hMin && STATUS(offset) == ACTIVE) {
      hq[hToIndex(hMin)].push(offset);
      break;
    }

  node_at_level[hToIndex(hMin)] = true;

  this->flood(hToIndex(hMin));

  Node* root = index[hToIndex(hMin)][0];

  // crop STATUS image to recover original dimensions

  this->m_parent->STATUS =
      this->STATUS.crop(back[0], this->STATUS.getSizeX() - front[0], back[1],
                        this->STATUS.getSizeY() - front[1], back[2],
                        this->STATUS.getSizeZ() - front[2]);

  this->m_parent->index = this->index;
  this->m_parent->hMin = this->hMin;

  return root;
}

// initialize global index for nodes
template <class T>
void SalembierRecursiveImplementation<T>::init(Image<T>& img,
                                               FlatSE& connexity) {
  FlatSE se = connexity;

  const TSize* tmpSize = img.getSize();
  const TCoord* tmpBack = se.getNegativeOffsets();
  const TCoord* tmpFront = se.getPositiveOffsets();

  for (int i = 0; i <= 2; i++) {
    oriSize[i] = tmpSize[i];
    back[i] = tmpBack[i];
    front[i] = tmpFront[i];
  }

  imBorder = img;
  imGradient = morphologicalGradient(img, connexity);
  STATUS.setSize(img.getSize());
  STATUS.fill(ACTIVE);

  addBorders(imBorder, back, front, BORDER);
  addBorders(STATUS, back, front, BORDER_STATUS);
  se.setContext(imBorder.getSize());

  this->hMin = img.getMin();
  this->hMax = img.getMax();
  this->numberOfLevels = hMax - hMin + 1;

  index.resize(numberOfLevels);

  hq = new std::queue<TOffset>[numberOfLevels];

  // we take a (max-min+1) * (number of grey-levels at level h)
  //  so we compute histogram

  int* histo;
  histo = new int[numberOfLevels];
  for (int i = 0; i < numberOfLevels; i++) histo[i] = 0;

  typename Image<T>::iterator it;
  typename Image<T>::iterator end = img.end();

  for (it = img.begin(); it != end; ++it) {
    int n = *it - hMin;
    assert(n >= 0 && n < numberOfLevels);
    histo[n]++;
  }

  for (int i = 0; i < numberOfLevels; i++) {
    int sizeOfLevelH = histo[i];

    index[i].resize(sizeOfLevelH);

    for (int j = 0; j < sizeOfLevelH; j++) index[i][j] = 0;
  }

  this->se = se;

  this->number_nodes.resize(numberOfLevels);
  this->node_at_level.resize(numberOfLevels);

  for (int i = 0; i < numberOfLevels; i++) {
    this->number_nodes[i] = 0;
    this->node_at_level[i] = false;
  }

  delete[] histo;
}

template <class T>
void SalembierRecursiveImplementation<T>::link_node(Node* tree, Node* child) {
  child->father = tree;
  tree->childs.push_back(child);
}

template <class T>
Node* SalembierRecursiveImplementation<T>::new_node(int h, int n) {
  Node* res = new Node;
  res->ori_h = h;
  res->h = h;
  res->label = n;

  return res;
}

}  // namespace LibTIM
