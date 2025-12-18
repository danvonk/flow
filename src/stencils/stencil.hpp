#pragma once

#include "config.hpp"
#include "fields/flow_field.hpp"

#include <concepts>

template <typename Stencil>
concept FieldStencil2D =
    requires(Stencil s, FlowFieldView &field, int i, int j) {
      { s(field, i, j) } -> std::same_as<void>;
    };

template <typename Stencil>
concept BoundaryStencil2D =
    requires(Stencil s, FlowFieldView &field, int i, int j) {
      { s.apply_left_wall(field, i, j) } -> std::same_as<void>;
      { s.apply_right_wall(field, i, j) } -> std::same_as<void>;
      { s.apply_top_wall(field, i, j) } -> std::same_as<void>;
      { s.apply_bottom_wall(field, i, j) } -> std::same_as<void>;
    };
