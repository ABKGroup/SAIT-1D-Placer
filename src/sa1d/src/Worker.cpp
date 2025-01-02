/////////////////////////////////////////////////////////////////////////////
// Authors: Sayak Kundu (sakundu@ucsd.edu), Zhiang Wang (zhw033@ucsd.edu)
//          Dooseok Yoon (d3yoon@ucsd.edu)
// Copyright (c) 2024, The Regents of the University of California
// All rights reserved.
//
// BSD 3-Clause License
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////

#include "Worker.h"
#include <cstdint>
#include <limits>
#include <numeric>
#include "odb/dbTypes.h"
#include "sa1d/Objects.h"
#include "sa1d/OptSA.h"
#include <iostream>


namespace sa1d {

using utl::SA1D;

SAWorker::SAWorker(OptSA* opt, Logger* logger, int worker_id)
  : opt_(opt),  logger_(logger), worker_id_(worker_id),
    cells_(opt->getCells()), nets_(opt->getNets())
{
  // Constructor implementation
  // Initialize other member variables if needed
  total_cells_ = static_cast<int>(cells_.size());
  num_nets_ = static_cast<int>(nets_.size());
  net_hpwl_.resize(num_nets_);
  affected_cells_.reserve(total_cells_);
}

SAWorker::~SAWorker()
{
  // Destructor implementation
  // Clean up resources if necessary
}

void SAWorker::setRandomSeed(int seed) 
{
  seed_ = seed;
  rng_.seed(seed_);
}

void SAWorker::setTemp(float temp)
{
  temp_ = temp;
}

void SAWorker::setCoolingRate(float cooling_rate)
{
  cooling_rate_ = cooling_rate;
}

void SAWorker::setNumMovePerIter(int num_move_per_iter)
{
  num_move_per_iter_ = num_move_per_iter;
}

void SAWorker::setMaxIter(int max_iter)
{
  max_iter_ = max_iter;
}

void SAWorker::setMoveProbs(const std::vector<float>& move_probs)
{
  action_prob_ = move_probs;
}

void SAWorker::initCellOrderRandom()
{
  cell_order_.resize(total_cells_);
  std::iota(cell_order_.begin(), cell_order_.end(), 0);
  std::shuffle(cell_order_.begin(), cell_order_.end(), rng_);
  std::vector<odb::dbOrientType> cell_orient;
  cell_orient.reserve(total_cells_);
  auto default_orient = odb::dbOrientType::R0;
  for (int i = 0; i < total_cells_; i++) {
    if (distribution_(rng_) <= 0.5) {
      cell_orient.push_back(default_orient);
    } else {
      cell_orient.push_back(orientMirrorY(default_orient));
    }
  }
  updateCellLocs(cell_orient);
  norm_hpwl_ = static_cast<float>(getTotalHPWL());
  logger_->report("[INFO] Norm HPWL = {}", static_cast<int>(norm_hpwl_));
}


void SAWorker::initCellOrder(
  const std::vector<int>& cell_order, 
  const std::vector<odb::dbOrientType>& cell_orients)
{
  cell_order_ = cell_order;
  updateCellLocs(cell_orients);
  norm_hpwl_ = static_cast<float>(getTotalHPWL());
  logger_->report("[INFO] Norm HPWL = {}", static_cast<int>(norm_hpwl_));
}

void SAWorker::updateCellLocs(const std::vector<odb::dbOrientType>& orients)
{
  cell_locs_.clear();
  cell_locs_.resize(total_cells_);
  int prev_width = opt_->getCoreLLX();
  int y = opt_->getCoreLLY();
  int x = 0.0;
  for (int i = 0; i < total_cells_; i++) {
    int cell_id = cell_order_[i];
    const Cell& cell = cells_[cell_id];
    x += prev_width;
    cell_locs_[cell_id] = CellLoc(x, y, orients[cell_id]);
    prev_width = cell.getWidth();
  }
  updateNetHPWL();
}

// We only consider x dim (y dim is fixed)
void SAWorker::calcNetHPWL(int net_id)
{
  int lx = std::numeric_limits<int>::max();
  int ux = std::numeric_limits<int>::min();  
  const auto& net = nets_[net_id];
  if (net.bterm_flag == true) {
    lx = std::min(lx, net.bterm_box.xMin());
    ux = std::max(ux, net.bterm_box.xMax());
  }  

  for (const auto& term : net.getTerms()) {
    int cell_id = term.cell_id;
    int mterm_id = term.mterm_id;
    const auto& cell_loc = cell_locs_[cell_id];
    auto pin_loc = cells_[cell_id].getPinLocation(
      cell_loc.x, cell_loc.y, mterm_id, cell_loc.orient);
    lx = std::min(lx, pin_loc.first);      
    ux = std::max(ux, pin_loc.first);
  }

  net_hpwl_[net_id] = netHPWL(ux - lx, lx, ux);
}

void SAWorker::updateNetHPWL()
{
  for (int net_id = 0; net_id < num_nets_; net_id++) {
    calcNetHPWL(net_id);
  } 
}

int64_t SAWorker::getTotalHPWL()
{
  updateNetHPWL();
  int64_t total_hpwl = 0;
  for (const auto& net : net_hpwl_) {
    total_hpwl += net.hpwl;
  }
  return total_hpwl;
}


void SAWorker::checkCellOrder()
{
  int prev_width = opt_->getCoreLLX();
  int expected_x = 0;
  for (int i = 0; i < total_cells_; i++) {
    int cell_id = cell_order_[i];
    const Cell& cell = cells_[cell_id];
    expected_x += prev_width;
    int x = cell_locs_[cell_id].x;
    if (x != expected_x) {
      logger_->error(sa1d::SA1D, 1, "Cell order mismatch at index: {} (cell id = {}) "
                      "Expected x: {} Found x: {}", i, cell_id, expected_x, x);
    }
    prev_width = cell.getWidth();
  }
}


int SAWorker::selectCell()
{
  int cell_id = static_cast<int>(distribution_(rng_) * total_cells_);
  return cell_id;
}


int SAWorker::selectMoveType()
{
  int move_type = static_cast<int>(distribution_(rng_) * move_count_);
  return move_type;  
}


/*
int SAWorker::swapCells(int id1, int id2)
{
  if (id1 == id2) {
    return 0;
  }

  // id1 is less than id2
  if (id1 > id2) {
    std::swap(id1, id2);
  }

  bool isSameWidth = opt_->cells_[cell_order_[id1]].getWidth() ==
                    opt_->cells_[cell_order_[id2]].getWidth();
  // Check if the cell width is the same or not
  CellLocMap moved_cells;
  // Calculate the HPWL of the affected nets
  int64_t update_hpwl = 0;
  int64_t current_hpwl = 0;
  set<Net*> affected_nets;
  if (isSameWidth) {
    for (auto net : opt_->cells_[cell_order_[id1]].nets_) {
      affected_nets.insert(net);
    }
    for (auto net : opt_->cells_[cell_order_[id2]].nets_) {
      affected_nets.insert(net);
    }

    moved_cells[&opt_->cells_[cell_order_[id1]]] = 
          cell_locs_[&opt_->cells_[cell_order_[id2]]];
    moved_cells[&opt_->cells_[cell_order_[id2]]] = 
          cell_locs_[&opt_->cells_[cell_order_[id1]]];
  } else {
    int width1 = opt_->cells_[cell_order_[id1]].getWidth();
    int width2 = opt_->cells_[cell_order_[id2]].getWidth();
    int delta_width = width1 - width2;
    moved_cells[&opt_->cells_[cell_order_[id2]]] = 
          cell_locs_[&opt_->cells_[cell_order_[id1]]];
    moved_cells[&opt_->cells_[cell_order_[id1]]] = 
          cell_locs_[&opt_->cells_[cell_order_[id2]]];
    moved_cells[&opt_->cells_[cell_order_[id1]]].first.first -= delta_width;
    for ( int i = id1+1; i < id2; i++ ) {
      moved_cells[&opt_->cells_[cell_order_[i]]] = 
            cell_locs_[&opt_->cells_[cell_order_[i]]];
      moved_cells[&opt_->cells_[cell_order_[i]]].first.first -= delta_width;
    }

    for ( int i = id1; i <= id2; i++ ) {
      for (auto net : opt_->cells_[cell_order_[i]].nets_) {
        affected_nets.insert(net);
      }
    }
  }

  for (auto net : affected_nets) {
    update_hpwl += net->getDeltaHPWL(cell_locs_, moved_cells);
    current_hpwl += net_hpwl_[net];
  }
  int delta_hpwl = static_cast<int>(update_hpwl - current_hpwl);
  return delta_hpwl;
}

void SAWorker::swapCellsAccept(int id1, int id2)
{
  if (id1 == id2) {
    return;
  }

  // id1 is less than id2
  if (id1 > id2) {
    std::swap(id1, id2);
  }

  bool isSameWidth = opt_->cells_[cell_order_[id1]].getWidth() ==
                    opt_->cells_[cell_order_[id2]].getWidth();
  
  // Calculate the HPWL of the affected nets
  set<Net*> affected_nets;
  if (isSameWidth) {
    for (auto net : opt_->cells_[cell_order_[id1]].nets_) {
      affected_nets.insert(net);
    }
    for (auto net : opt_->cells_[cell_order_[id2]].nets_) {
      affected_nets.insert(net);
    }
    std::swap(cell_order_[id1], cell_order_[id2]); 
    // Updated cell_locs_ based on the swap
    std::swap(cell_locs_[&opt_->cells_[cell_order_[id1]]], 
            cell_locs_[&opt_->cells_[cell_order_[id2]]]);
  } else {
    int width1 = opt_->cells_[cell_order_[id1]].getWidth();
    int width2 = opt_->cells_[cell_order_[id2]].getWidth();
    int delta_width = width1 - width2;
    std::swap(cell_order_[id1], cell_order_[id2]);
    std::swap(cell_locs_[&opt_->cells_[cell_order_[id1]]], 
            cell_locs_[&opt_->cells_[cell_order_[id2]]]);
    for ( int i = id1+1; i <= id2; i++ ) {
      cell_locs_[&opt_->cells_[cell_order_[i]]].first.first -= delta_width;
    }
    for ( int i = id1; i <= id2; i++ ) {
      for (auto net : opt_->cells_[cell_order_[i]].nets_) {
        affected_nets.insert(net);
      }
    }
  }
  // Update net HPWL of the affected nets
  for (auto net : affected_nets) {
    net_hpwl_[net] = net->getHPWL(cell_locs_);
  }
}

int SAWorker::moveCell(int id, int n)
{
  // If new location is same as the current location return 0
  if (n == id ) {
    return 0;
  }

  // Get the cell width
  int width1 = opt_->cells_[cell_order_[id]].getWidth();
  int width2 = opt_->cells_[cell_order_[n]].getWidth();
  
  int new_y = cell_locs_[&opt_->cells_[cell_order_[n]]].first.second;
  // Match the left edge of the cell
  int new_x = cell_locs_[&opt_->cells_[cell_order_[n]]].first.first;
  // Get the new location
  if ( n > id ) {
    // Match the right edge of the cell
    new_x = cell_locs_[&opt_->cells_[cell_order_[n]]].first.first + width2 - width1;
  }
  dbOrientType orient = cell_locs_[&opt_->cells_[cell_order_[id]]].second;

  CellLocMap moved_cells;
  set<Net*> affected_nets;
  
  // Update location of the moved cell
  moved_cells[&opt_->cells_[cell_order_[id]]] = {{new_x, new_y}, orient};
  for (auto net : opt_->cells_[cell_order_[id]].nets_) {
    affected_nets.insert(net);
  }
  
  // Update location of the cells in between
  int id1, id2, direction;
  int width = width1;
  if ( id < n ) {
    id1 = id + 1;
    id2 = n;
    direction = -1;
  } else {
    id1 = n;
    id2 = id - 1;
    direction = 1;
  }

  for ( int i = id1; i <= id2; i++ ) {
    moved_cells[&opt_->cells_[cell_order_[i]]] = 
          cell_locs_[&opt_->cells_[cell_order_[i]]];
    moved_cells[&opt_->cells_[cell_order_[i]]].first.first += direction * width;
    for (auto net : opt_->cells_[cell_order_[i]].nets_) {
      affected_nets.insert(net);
    }
  }

  // Calculate the HPWL of the affected nets
  int64_t update_hpwl = 0;
  int64_t current_hpwl = 0;
  for (auto net : affected_nets) {
    update_hpwl += net->getDeltaHPWL(cell_locs_, moved_cells);
    current_hpwl += net_hpwl_[net];
  }
  int delta_hpwl = static_cast<int>(update_hpwl - current_hpwl);
  return delta_hpwl;
}

// Move the cell at index id to index n
void SAWorker::moveCellAccept(int id, int n)
{
  // If new location is the same as the current location, return
  if (n == id) {
    return;
  }

  // Identify the cell being moved
  int cell_id = cell_order_[id];
  int width_src = opt_->cells_[cell_id].getWidth();
  int width_dest = opt_->cells_[cell_order_[n]].getWidth();

  // Compute new location for the moved cell
  int new_x = cell_locs_[&opt_->cells_[cell_order_[n]]].first.first;
  int new_y = cell_locs_[&opt_->cells_[cell_order_[n]]].first.second;
  auto orient = cell_locs_[&opt_->cells_[cell_id]].second;

  // Collect affected nets (from the cell that’s moving)
  set<Net*> affected_nets;
  for (auto net : opt_->cells_[cell_id].nets_) {
    affected_nets.insert(net);
  }

  // Shift intermediate cells
  if (id < n) {
    // Shift cells in [id+1..n] one position left
    for (int i = id; i < n; i++) {
      cell_order_[i] = cell_order_[i + 1];
      int shifted_cell_id = cell_order_[i];
      
      // Update location of each shifted cell
      cell_locs_[&opt_->cells_[shifted_cell_id]].first.first -= width_src;
      
      // Mark nets as affected
      for (auto net : opt_->cells_[shifted_cell_id].nets_) {
        affected_nets.insert(net);
      }
    }
  } else {
    // Shift cells in [n..id-1] one position right
    for (int i = id; i > n; i--) {
      cell_order_[i] = cell_order_[i - 1];
      int shifted_cell_id = cell_order_[i];
      // Update location of each shifted cell
      cell_locs_[&opt_->cells_[shifted_cell_id]].first.first += width_src;
      // Mark nets as affected
      for (auto net : opt_->cells_[shifted_cell_id].nets_) {
        affected_nets.insert(net);
      }
    }
  }

  // Place the moved cell at position n
  cell_order_[n] = cell_id;

  // If moving to the right, align the right edge
  // (destination width might be larger or smaller than source width)
  if (n > id) {
    // Align right edges: new_x = dest's left X + (destWidth - srcWidth)
    new_x += (width_dest - width_src);
  }

  // Update location of the moved cell
  cell_locs_[&opt_->cells_[cell_id]] = {{new_x, new_y}, orient};

  // Update net HPWL of all affected nets
  for (auto net : affected_nets) {
    net_hpwl_[net] = net->getHPWL(cell_locs_);
  }
}
*/

int SAWorker::swapCells(int id1, int id2)
{
  if (id1 == id2) {
    return 0;
  }

  // id1 is less than id2
  if (id1 > id2) {
    std::swap(id1, id2);
  }

  int cell_id1 = cell_order_[id1];
  int cell_id2 = cell_order_[id2];

  auto& cell1 = cells_[cell_order_[id1]];
  auto& cell2 = cells_[cell_order_[id2]];
  bool isSameWidth = cell1.getWidth() == cell2.getWidth();

  affected_cells_.clear();
  affected_nets_.clear();
  cell_locs_[cell_id1].backupX();
  cell_locs_[cell_id2].backupX();
  std::swap(cell_locs_[cell_id1].x, cell_locs_[cell_id2].x);
  std::swap(cell_order_[id1], cell_order_[id2]);
  affected_cells_.push_back(cell_id1);
  affected_cells_.push_back(cell_id2);
  
  if (isSameWidth == false) {
    int delta_width = cell1.getWidth() - cell2.getWidth();
    cell_locs_[cell_id1].x -= delta_width;
    for (int i = id1+1; i < id2; i++ ) {
      cell_locs_[cell_order_[i]].backupX();
      cell_locs_[cell_order_[i]].x -= delta_width;
      affected_cells_.push_back(cell_order_[i]);
    }
  }

  // update the affected nets
  for (auto& cell_id : affected_cells_) {
    for (auto& net_id : cells_[cell_id].nets) {
      if (net_id == -1) {
        continue;
      }
      affected_nets_.insert(net_id);
    }
  }

  for (auto& net_id : affected_nets_) {
    net_hpwl_[net_id].backup();
  }

  for (auto& cell_id : affected_cells_) {
    updateDeltaHPWLUtil(cell_id);
  }

  int delta_hpwl = 0;
  for (auto& net_id : affected_nets_) {
    delta_hpwl += updateDeltaHPWL(net_id);
  }

  return delta_hpwl;
}


void SAWorker::swapCellsReject(int id1, int id2)
{
  for (auto& cell_id : affected_cells_) {
    cell_locs_[cell_id].restoreX();
  }
  
  for (auto& net_id : affected_nets_) {
    net_hpwl_[net_id].restore();
  }
  
  std::swap(cell_order_[id1], cell_order_[id2]);
}


int SAWorker::moveCell(int id, int n)
{
  // If new location is same as the current location return 0
  if (n == id ) {
    return 0;
  }

  // Get the cell width
  int cell_id1 = cell_order_[id];
  int cell_id2 = cell_order_[n];

  int width1 = cells_[cell_id1].getWidth();
  int width2 = cells_[cell_id2].getWidth();
  cell_locs_[cell_id1].backupX();
  cell_locs_[cell_id2].backupX();

  affected_cells_.clear();
  affected_nets_.clear();
  // We need to consider two possible cases:
  // if id > n, we need to insert the cell in front of n and shift the cells in between towards the right
  // if id < n, we need to insert the cell behind n and shift the cells in between towards the left
  if (id > n) {
    int new_x = cell_locs_[cell_id2].x; // Match the left edge of the cell
    cell_locs_[cell_id1].x = new_x;
    for (int i = n; i < id; i++) {
      auto& cell_id_temp = cell_order_[i];
      cell_locs_[cell_id_temp].backupX();
      cell_locs_[cell_id_temp].x += width1;
    }
    for (int i = id; i > n; i--) {
      cell_order_[i] = cell_order_[i - 1];
      affected_cells_.push_back(cell_order_[i]);
    }
    cell_order_[n] = cell_id1;
    affected_cells_.push_back(cell_id1);
  } else { // id < n
    int new_x = cell_locs_[cell_id2].x + width2 - width1; // Match the right edge of the cell
    cell_locs_[cell_id1].x = new_x;
    for (int i = id + 1; i <= n; i++) {
      auto& cell_id_temp = cell_order_[i];
      cell_locs_[cell_id_temp].backupX();
      cell_locs_[cell_id_temp].x -= width1;
    }
    for (int i = id; i < n; i++) {
      cell_order_[i] = cell_order_[i + 1];
      affected_cells_.push_back(cell_order_[i]);
    }
    cell_order_[n] = cell_id1;
    affected_cells_.push_back(cell_id1);
  }

  // update the affected nets
  // update the affected nets
  for (auto& cell_id : affected_cells_) {
    for (auto& net_id : cells_[cell_id].nets) {
      if (net_id == -1) {
        continue;
      }
      affected_nets_.insert(net_id);
    }
  }

  for (auto& net_id : affected_nets_) {
    net_hpwl_[net_id].backup();
  }

  for (auto& cell_id : affected_cells_) {
    updateDeltaHPWLUtil(cell_id);
  }

  int delta_hpwl = 0;
  for (auto& net_id : affected_nets_) {
    delta_hpwl += updateDeltaHPWL(net_id);
  }

  return delta_hpwl;
}


void SAWorker::moveCellReject(int id, int n)
{
  // If new location is same as the current location return 0
  if (n == id ) {
    return;
  }
  
  for (auto& cell_id : affected_cells_) {
    cell_locs_[cell_id].restoreX();
  }
  
  for (auto& net_id : affected_nets_) {
    net_hpwl_[net_id].restore();
  }
   
  int cell_id2 = cell_order_[n];
  if (id > n) {
    for (int i = n; i < id; i++) {
      cell_order_[i] = cell_order_[i + 1];
    }
    cell_order_[id] = cell_id2;    
  } else { // id < n
    for (int i = n; i > id; i--) {
      cell_order_[i] = cell_order_[i - 1];
    }
    cell_order_[id] = cell_id2;
  } 
}


// Key point: use the boundary condition to calculate the delta HPWL
void SAWorker::updateDeltaHPWLUtil(int cell_id)
{
  auto& cell_loc = cell_locs_[cell_id];
  const auto& cell_x = cell_loc.x;
  const auto& cell_y = cell_loc.y;  
  const auto& orient = cell_loc.orient;

  const auto& cell_pre_x = cell_locs_[cell_id].pre_x;
  auto& cell = cells_[cell_id];
  for (int mterm_id = 0; mterm_id < cell.getMtermCount(); mterm_id++) {
    auto mterm_loc = cell.getPinLocation(cell_x, cell_y, mterm_id, orient);
    auto pre_mterm_loc = cell.getPinLocation(cell_pre_x, cell_y, mterm_id, orient);
    int net_id = cell.nets[mterm_id];
    if (net_id == -1)  {
      continue;
    }

    auto& net = net_hpwl_[net_id];
    if (net.scratch_flag == true) {
      continue;
    }  

    if (mterm_loc.first <= net.lx && pre_mterm_loc.first < net.ux) {
      net.lx = mterm_loc.first;
    } else if (mterm_loc.first >= net.ux && pre_mterm_loc.first > net.lx) {
      net.ux = mterm_loc.first;
    } else {
      net.scratch_flag = true;  
    }
  }
}

// We only consider x dim (y dim is fixed)
int SAWorker::updateDeltaHPWL(int net_id)
{
  auto& net_hpwl = net_hpwl_[net_id];
  if (net_hpwl.scratch_flag == false) {
    net_hpwl.hpwl = net_hpwl.ux - net_hpwl.lx;
  } else {
    int lx = std::numeric_limits<int>::max();
    int ux = std::numeric_limits<int>::min();  
    const auto& net = nets_[net_id];
    if (net.bterm_flag == true) {
      lx = std::min(lx, net.bterm_box.xMin());
      ux = std::max(ux, net.bterm_box.xMax());
    }  
    for (const auto& term : net.getTerms()) {
      int cell_id = term.cell_id;
      int mterm_id = term.mterm_id;
      const auto& cell_loc = cell_locs_[cell_id];
      auto pin_loc = cells_[cell_id].getPinLocation(
        cell_loc.x, cell_loc.y, mterm_id, cell_loc.orient);
      lx = std::min(lx, pin_loc.first);      
      ux = std::max(ux, pin_loc.first);
    }
    net_hpwl.ux = ux;
    net_hpwl.lx = lx;
    net_hpwl.hpwl = ux - lx;
  }
  return net_hpwl.getDeltaHPWL();
}


int SAWorker::flipCell(int id)
{
  int cell_id = cell_order_[id];
  cell_locs_[cell_id].backupOrient();
  cell_locs_[cell_id].flipOrient();
  affected_nets_.clear();

  // To speed up the calculation of delta HPWL,
  // we calculate the change of HPWL in a cell-based manner
  for (auto net_id : cells_[cell_id].getNets()) {
    if (net_id == -1) {
      continue;
    }
    
    affected_nets_.insert(net_id);
    net_hpwl_[net_id].backup();
    updateDeltaHPWLUtil(cell_id);
  }

  int delta_hpwl = 0;
  for (auto& net_id : affected_nets_) {
    delta_hpwl += updateDeltaHPWL(net_id);
  }

  return delta_hpwl;
}


void SAWorker::flipCellReject(int id)
{
  int cell_id = cell_order_[id];
  cell_locs_[cell_id].restoreOrient();
  for (auto net_id : affected_nets_) {
    net_hpwl_[net_id].restore();
  }
}


// define the same guardband
float safeExp(float x) 
{
  x = std::min(x, 50.0f);
  x = std::max(x, -50.0f);
  return std::exp(x);  
}

// Define the guard for exp to avoid overflow
void SAWorker::run()
{
  accept_count_ = 0;
  int64_t curr_hpwl = getTotalHPWL();
  for (int iter = 0; iter < max_iter_; iter++) {
    for (int move_id = 0; move_id < num_move_per_iter_; move_id++) {
      int move_type = selectMoveType();
      int id1 = selectCell();
      int id2 = selectCell();
      int delta_hpwl = 0;
      switch (move_type) {
        case 0:
          delta_hpwl = swapCells(id1, id2);
          break;
        case 1:
          delta_hpwl = moveCell(id1, id2);
          break;
        case 2:
          delta_hpwl = flipCell(id1);
          break;
        default:
          logger_->error(sa1d::SA1D, 4, "Invalid move type selected.");
      }

      bool accept_flag = delta_hpwl <= 0.0 || 
        distribution_(rng_) < safeExp(-1.0 * delta_hpwl / norm_hpwl_ / temp_);
      // bool accept_flag = true;

      if (accept_flag == false) {
        // Reject the move
        switch (move_type) {
          case 0:
            swapCellsReject(id1, id2);
            break;
          case 1:
            moveCellReject(id1, id2);
            break;
          case 2:
            flipCellReject(id1);
            break;
          default:
            break;
        }
      } else {
        accept_count_++;
        curr_hpwl += delta_hpwl;
      }

      if (false) {// for debug only
        auto net_hpwl = net_hpwl_; 
        if (curr_hpwl != getTotalHPWL()) {
          std::cout << "iter = " << iter << " move_id = " << move_id << " selected cell id = " << id1 << " " << id2 << " "
                    << "accept_flag = " << accept_flag << std::endl;
          for (int i = 0; i < num_nets_; i++) {
            if (net_hpwl_[i].hpwl != net_hpwl[i].hpwl) {
              std::cout << "net hpwl mismatch at index = " << i << " " << net_hpwl_[i].hpwl << " " << net_hpwl[i].hpwl << " "
                        << "new lx = " << net_hpwl_[i].lx << " old lx =  " << net_hpwl[i].lx << " "
                        << "new ux = " << net_hpwl_[i].ux << " old ux =  " << net_hpwl[i].ux << std::endl;
            }
          }
                    
          logger_->error(sa1d::SA1D, 5, "curr_hpwl does not match !!!");
        }

        auto cell_order = getFinalOrdering();
        auto orients = getFinalOrients();
        initCellOrder(cell_order, orients);
        if (curr_hpwl != getTotalHPWL()) {
          std::cout << "iter = " << iter << " move_id = " << move_id << " selected cell id = " << id1 << " " << id2 << " "
                    << "accept_flag = " << accept_flag << std::endl;
          logger_->error(sa1d::SA1D, 6, "the updated location or flip operation is wrong !!!");
        }

          /*
          // for testing
          const std::vector<int>& cell_order = getFinalOrdering();
          const std::vector<odb::dbOrientType>& orients = getFinalOrients();
          auto net_hpwl = net_hpwl_;      
          auto cell_loc = cell_locs_;

          initCellOrder(cell_order, orients);
          if (curr_hpwl != getTotalHPWL()) {
            auto new_orient = getFinalOrients();
            for (int i = 0; i < total_cells_; i++) {
              if (orients[i] != new_orient[i]) {
                std::cout << "orient mismatch at index = " << i << " " << orients[i] << " " << new_orient[i] << std::endl;
              }

              if (cell_order_[i] != cell_order[i]) {
                std::cout << "cell order mismatch at index = " << i << " " << cell_order_[i] << " " << cell_order[i] << std::endl;
              }
            }

            for (int i = 0; i < total_cells_; i++) {
              if (cell_loc[i].x != cell_locs_[i].x) {
                std::cout << "cell x mismatch at index = " << i << " " << cell_loc[i].x << " " << cell_locs_[i].x << std::endl;
              }

              if (cell_loc[i].orient != cell_locs_[i].orient) {
                std::cout << "cell orient mismatch at index = " << i << " " << cell_loc[i].orient << " " << cell_locs_[i].orient << std::endl;
              }
            }

            for (int i = 0; i < num_nets_; i++) {
              if (net_hpwl_[i].hpwl != net_hpwl[i].hpwl) {
                std::cout << "net hpwl mismatch at index = " << i << " " << net_hpwl_[i].hpwl << " " << net_hpwl[i].hpwl << std::endl;
              }
            }
            
            std::cout << "iter = " <<  iter << " move_id = " << move_id << " selected cell id = " << id1 << " " << id2 << std::endl;
         
          }
          */
      }
    }
    temp_ = temp_ * cooling_rate_;
  }
}


void SAWorker::reportDetails() {
  /*
  // Report Current Iteration and Temperature
  logger_->report("###########################################");
  logger_->report("## Worker ID: {:<26} ##", worker_id_);
  logger_->report("## Current Iteration: {:<18} ##", iter_);
  
  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> elapsed_time = 
      end_time - iter_start_time_;
  // Duration in seconds up to 3 decimal places
  double elapsed_time_sec = elapsed_time.count() / 1000.0;
  logger_->report("## Elapsed Time for Iteration: {:<9.3f} ##", 
                  elapsed_time_sec);

  // Report initialization time iter_start_time_ - start_time_
  elapsed_time = iter_start_time_ - start_time_;
  double init_time_sec = elapsed_time.count() / 1000.0;
  logger_->report("## Initialization Time: {:<16.3f} ##", init_time_sec);

  // Report Move Probs, Max Temps, Min Temp, Max Iter, Seed
  logger_->report("## Move Probs: {:<4.2f} {:<4.2f} {:<15.2f} ##", 
                  action_prob_[0], action_prob_[1], action_prob_[2]);
  logger_->report("## Accepted Moves: {:<21} ##", accept_count_);
  logger_->report("## Max Temps: {:<26.2f} ##", max_temp_);
  logger_->report("## Current Temp: {:<23.2f} ##", temp_);
  logger_->report("## Min Temp: {:<27.2f} ##", min_temp_);
  logger_->report("## Max Iter: {:<27} ##", max_iter_);
  logger_->report("## Seed: {:<31} ##", seed_);
  
  // updateNetHPWL();
  int64_t total_hpwl = getTotalHPWL();

  logger_->report("## Total HPWL: {:<25.3f} ##",
                 opt_->block_->dbuToMicrons(total_hpwl));
  logger_->report("###########################################");
  */
  checkCellOrder();
}

}  // namespace sa1d
