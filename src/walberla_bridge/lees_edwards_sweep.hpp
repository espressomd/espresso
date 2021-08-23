#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Environment.h"
#include "core/timing/RemainingTimeLogger.h"

#include "field/AddToStorage.h"
#include "field/vtk/VTKWriter.h"

#include "timeloop/SweepTimeloop.h"

#include <memory>
#include <cmath>
/*
using namespace walberla;

using flag_t      = walberla::uint8_t;
using FlagField_T = FlagField< flag_t >;

class LeesEdwardsUpdate
{
 public:
   LeesEdwardsUpdate(const std::shared_ptr< StructuredBlockForest >& blocks, BlockDataID fieldID, real_t offset)
      : blocks_(blocks), fieldID_(fieldID), offset_(offset)
   {}

   void operator()(IBlock* block)
   {
      // TODO should dimension_x contain the ghost layers or not. At the moment value is 64 with GL it is 66. In the
      // lbmpy Leed Edwards this is mixed. Probably not good


      // Top cells
      if (blocks_->atDomainYMaxBorder(*block))
      {
         uint_t dimension_x = blocks_->getNumberOfXCells(*block);
         real_t weight      = fmod(offset_ + real_c(dimension_x), 1.0);

         auto pdf_field = block->getData< PdfField_T >(fieldID_);
         // auto pdf_field = (*bc).block->template getData<PdfField>(m_pdf_field_id);

         CellInterval ci;
         pdf_field->getGhostRegion(stencil::N, ci, 1, true);

         for (auto cell = ci.begin(); cell != ci.end(); ++cell)
         {
            cell_idx_t x = cell->x();

            uint_t ind1 = uint_c(floor(x - offset_)) % dimension_x;
            uint_t ind2 = uint_c(ceil(x - offset_)) % dimension_x;

            for (uint_t q = 0; q < Stencil_T::Q; ++q)
            {
               pdf_field->get(*cell, 0) = (1 - weight) * pdf_field->get(cell_idx_c(ind1), cell->y(), cell->z(), q) +
                                          weight * pdf_field->get(cell_idx_c(ind2), cell->y(), cell->z(), q);
            }
         }
      }
      // Bottom cells
      if (blocks_->atDomainYMinBorder(*block))
      {
         uint_t dimension_x = blocks_->getNumberOfXCells(*block);
         real_t weight      = fmod(offset_ + real_c(dimension_x), 1.0);

         auto pdf_field = block->getData< PdfField_T >(fieldID_);

         CellInterval ci;
         pdf_field->getGhostRegion(stencil::S, ci, 1, true);

         for (auto cell = ci.begin(); cell != ci.end(); ++cell)
         {
            cell_idx_t x = cell->x();

            uint_t ind1 = uint_c(floor(x + offset_)) % dimension_x;
            uint_t ind2 = uint_c(ceil(x + offset_)) % dimension_x;

            for (uint_t q = 0; q < Stencil_T::Q; ++q)
            {
               pdf_field->get(*cell, 0) = (1 - weight) * pdf_field->get(cell_idx_c(ind1), cell->y(), cell->z(), q) +
                                          weight * pdf_field->get(cell_idx_c(ind2), cell->y(), cell->z(), q);
            }
         }
      }
   }

 private:
   const shared_ptr< StructuredBlockForest >& blocks_;
   BlockDataID fieldID_;
   real_t offset_;
};

*/