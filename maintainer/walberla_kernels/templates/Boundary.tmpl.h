/*
 * Copyright (C) 2022-2023 The ESPResSo project
 * Copyright (C) 2020-2023 The waLBerla project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * Boundary class.
 * Adapted from the waLBerla source file
 * https://i10git.cs.fau.de/walberla/walberla/-/blob/fb076cd18daa6e2f24448349d1fffb974c845269/python/pystencils_walberla/templates/Boundary.tmpl.h
 */

#pragma once

#include <core/DataTypes.h>

{% if target is equalto 'cpu' -%}
#include <field/GhostLayerField.h>
{%- elif target is equalto 'gpu' -%}
#include <cuda/GPUField.h>
#include <cuda/FieldCopy.h>
{%- endif %}
#include <domain_decomposition/BlockDataID.h>
#include <domain_decomposition/IBlock.h>
#include <blockforest/StructuredBlockForest.h>
#include <field/FlagField.h>
#include <core/debug/Debug.h>

#include <functional>
#include <set>
#include <vector>

{% for header in interface_spec.headers %}
#include {{header}}
{% endfor %}

#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

namespace walberla {
namespace {{namespace}} {


class {{class_name}}
{
public:
    {{StructDeclaration|indent(4)}}


    class IndexVectors
    {
    public:
        using CpuIndexVector = std::vector<{{StructName}}>;

        enum Type {
            ALL = 0,
            INNER = 1,
            OUTER = 2,
            NUM_TYPES = 3
        };

        IndexVectors() = default;
        bool operator==(IndexVectors const &other) const { return other.cpuVectors_ == cpuVectors_; }

        {% if target == 'gpu' -%}
        ~IndexVectors() {
            for( auto & gpuVec: gpuVectors_)
                cudaFree( gpuVec );
        }
        {% endif -%}

        CpuIndexVector & indexVector(Type t) { return cpuVectors_[t]; }
        {{StructName}} * pointerCpu(Type t)  { return cpuVectors_[t].data(); }

        {% if target == 'gpu' -%}
        {{StructName}} * pointerGpu(Type t)  { return gpuVectors_[t]; }
        {% endif -%}

        void syncGPU()
        {
            {% if target == 'gpu' -%}
            for( auto & gpuVec: gpuVectors_)
                cudaFree( gpuVec );
            gpuVectors_.resize( cpuVectors_.size() );

            WALBERLA_ASSERT_EQUAL(cpuVectors_.size(), NUM_TYPES);
            for(size_t i=0; i < cpuVectors_.size(); ++i )
            {
                auto & gpuVec = gpuVectors_[i];
                auto & cpuVec = cpuVectors_[i];
                cudaMalloc( &gpuVec, sizeof({{StructName}}) * cpuVec.size() );
                cudaMemcpy( gpuVec, &cpuVec[0], sizeof({{StructName}}) * cpuVec.size(), cudaMemcpyHostToDevice );
            }
            {%- endif %}
        }

    private:
        std::vector<CpuIndexVector> cpuVectors_{NUM_TYPES};

        {% if target == 'gpu' -%}
        using GpuIndexVector = {{StructName}} *;
        std::vector<GpuIndexVector> gpuVectors_;
        {%- endif %}
    };

    {{class_name}}( const shared_ptr<StructuredBlockForest> & blocks,
                   {{kernel|generate_constructor_parameters(['indexVector', 'indexVectorSize'])}}{{additional_data_handler.constructor_arguments}})
        :{{additional_data_handler.initialiser_list}} {{ kernel|generate_constructor_initializer_list(['indexVector', 'indexVectorSize']) }}
    {
        auto createIdxVector = []( IBlock * const , StructuredBlockStorage * const ) { return new IndexVectors(); };
        indexVectorID = blocks->addStructuredBlockData< IndexVectors >( createIdxVector, "IndexField_{{class_name}}");
    };

    {{class_name}}({{kernel|generate_constructor_parameters(['indexVectorSize'])}}{{additional_data_handler.constructor_arguments}})
        : {{additional_data_handler.initialiser_list}} {{ kernel|generate_constructor_initializer_list(['indexVectorSize']) }}
    {};

    void run (
        {{- ["IBlock * block", kernel.kernel_selection_parameters, ["cudaStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}}
    );

    {% if generate_functor -%}
    void operator() (
        {{- ["IBlock * block", kernel.kernel_selection_parameters, ["cudaStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}}
    )
    {
        run( {{- ["block", kernel.kernel_selection_parameters, ["stream"] if target == 'gpu' else []] | identifier_list -}} );
    }
    {%- endif %}

    void inner (
        {{- ["IBlock * block", kernel.kernel_selection_parameters, ["cudaStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}}
    );

    void outer (
        {{- ["IBlock * block", kernel.kernel_selection_parameters, ["cudaStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}}
    );

    std::function<void (IBlock *)> getSweep( {{- [interface_spec.high_level_args, ["cudaStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}} )
    {
        return [ {{- ["this", interface_spec.high_level_args, ["stream"] if target == 'gpu' else []] | identifier_list -}} ]
               (IBlock * b)
               { this->run( {{- [ ['b'], interface_spec.mapping_codes, ["stream"] if target == 'gpu' else [] ] | identifier_list -}} ); };
    }

    std::function<void (IBlock *)> getInnerSweep( {{- [interface_spec.high_level_args, ["cudaStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}} )
    {
        return [ {{- [ ['this'], interface_spec.high_level_args, ["stream"] if target == 'gpu' else [] ] | identifier_list -}} ]
               (IBlock * b)
               { this->inner( {{- [ ['b'], interface_spec.mapping_codes, ["stream"] if target == 'gpu' else [] ] | identifier_list -}} ); };
    }

    std::function<void (IBlock *)> getOuterSweep( {{- [interface_spec.high_level_args, ["cudaStream_t stream = nullptr"] if target == 'gpu' else []] | type_identifier_list -}} )
    {
        return [ {{- [ ['this'], interface_spec.high_level_args, ["stream"] if target == 'gpu' else [] ] | identifier_list -}} ]
               (IBlock * b)
               { this->outer( {{- [ ['b'], interface_spec.mapping_codes, ["stream"] if target == 'gpu' else [] ] | identifier_list -}} ); };
    }

    template<typename FlagField_T>
    void fillFromFlagField( const shared_ptr<StructuredBlockForest> & blocks, ConstBlockDataID flagFieldID,
                            FlagUID boundaryFlagUID, FlagUID domainFlagUID)
    {
        for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
            fillFromFlagField<FlagField_T>({{additional_data_handler.additional_arguments_for_fill_function}}&*blockIt, flagFieldID, boundaryFlagUID, domainFlagUID );
    }


    template<typename FlagField_T>
    void fillFromFlagField({{additional_data_handler.additional_parameters_for_fill_function}}IBlock * block, ConstBlockDataID flagFieldID,
                            FlagUID boundaryFlagUID, FlagUID domainFlagUID )
    {
        auto * indexVectors = block->getData< IndexVectors > ( indexVectorID );
        auto & indexVectorAll = indexVectors->indexVector(IndexVectors::ALL);
        auto & indexVectorInner = indexVectors->indexVector(IndexVectors::INNER);
        auto & indexVectorOuter = indexVectors->indexVector(IndexVectors::OUTER);

        auto * flagField = block->getData< FlagField_T > ( flagFieldID );
        {{additional_data_handler.additional_field_data|indent(4)}}

        if( !(flagField->flagExists(boundaryFlagUID) && flagField->flagExists(domainFlagUID) ))
            return;

        auto boundaryFlag = flagField->getFlag(boundaryFlagUID);
        auto domainFlag = flagField->getFlag(domainFlagUID);

        auto inner = flagField->xyzSize();
        inner.expand( cell_idx_t(-1) );

        indexVectorAll.clear();
        indexVectorInner.clear();
        indexVectorOuter.clear();

        {% if inner_or_boundary -%}
        for( auto it = flagField->begin(); it != flagField->end(); ++it )
        {
            if( ! isFlagSet(it, domainFlag) )
                continue;
            {%- for dirIdx, dirVec, offset in additional_data_handler.stencil_info %}
            if ( isFlagSet( it.neighbor({{offset}} {%if dim == 3%}, 0 {%endif %}), boundaryFlag ) )
            {
                auto element = {{StructName}}(it.x(), it.y(), {%if dim == 3%} it.z(), {%endif %} {{dirIdx}} );
                {{additional_data_handler.data_initialisation(dirIdx)|indent(16)}}
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            {% endfor %}
        }
        {%else%}
        auto flagWithGLayers = flagField->xyzSizeWithGhostLayer();
        {% if single_link %}
        {{dtype}} dot = 0.0; {{dtype}} maxn = 0.0;
        cell_idx_t calculated_idx = 0;
        cell_idx_t dx = 0; cell_idx_t dy = 0; {%if dim == 3%}  cell_idx_t dz = 0; {% endif %}
        cell_idx_t sum_x = 0; cell_idx_t sum_y = 0; {%if dim == 3%} cell_idx_t sum_z = 0; {%endif %}
        {% endif -%}
        for( auto it = flagField->beginWithGhostLayerXYZ(); it != flagField->end(); ++it )
        {
            {% if single_link -%}
            sum_x = 0; sum_y = 0; {%if dim == 3%} sum_z = 0; {%endif %}
            {% endif %}
            if( ! isFlagSet(it, boundaryFlag) )
                continue;
            {%- for dirIdx, dirVec, offset in additional_data_handler.stencil_info %}
            if ( flagWithGLayers.contains(it.x() + cell_idx_c({{dirVec[0]}}), it.y() + cell_idx_c({{dirVec[1]}}), it.z() + cell_idx_c({{dirVec[2]}})) && isFlagSet( it.neighbor({{offset}} {%if dim == 3%}, 0 {%endif %}), domainFlag ) )
            {
                {% if single_link -%}
                sum_x += cell_idx_c({{dirVec[0]}}); sum_y += cell_idx_c({{dirVec[1]}}); {%if dim == 3%} sum_z += cell_idx_c({{dirVec[2]}}); {%endif %}
                {% else %}
                auto element = {{StructName}}(it.x(), it.y(), {%if dim == 3%} it.z(), {%endif %} {{dirIdx}} );
                {{additional_data_handler.data_initialisation(dirIdx)|indent(16)}}
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
                {% endif %}
            }
            {% endfor %}

        {% if single_link %}
            dot = 0.0; maxn = 0.0; calculated_idx = 0;
            if(sum_x != 0 or sum_y !=0 {%if dim == 3%} or sum_z !=0 {%endif %})
            {
            {%- for dirIdx, dirVec, offset in additional_data_handler.stencil_info %}
                dx = {{dirVec[0]}}; dy = {{dirVec[1]}}; {%if dim == 3%} dz = {{dirVec[2]}}; {% endif %}
                dot = numeric_cast< {{dtype}} >( dx*sum_x + dy*sum_y {%if dim == 3%} + dz*sum_z {% endif %});
                if (dot > maxn)
                {
                    maxn = dot;
                    calculated_idx = {{dirIdx}};
                }
            {% endfor %}
                auto element = {{StructName}}(it.x(), it.y(), {%if dim == 3%} it.z(), {%endif %} calculated_idx );
                {{additional_data_handler.data_initialisation(dirIdx)|indent(16)}}
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                indexVectorInner.push_back( element );
                else
                indexVectorOuter.push_back( element );
            }
        {% endif -%}

        }
        {% endif %}

        indexVectors->syncGPU();
    }

private:
    void run_impl(
        {{- ["IBlock * block", "IndexVectors::Type type",
             kernel.kernel_selection_parameters, ["cudaStream_t stream = nullptr"] if target == 'gpu' else []]
            | type_identifier_list -}}
   );

    BlockDataID indexVectorID;
    {{additional_data_handler.additional_member_variable|indent(4)}}
public:
    {{kernel|generate_members(('indexVector', 'indexVectorSize'))|indent(4)}}
};

} // namespace {{namespace}}
} // namespace walberla
