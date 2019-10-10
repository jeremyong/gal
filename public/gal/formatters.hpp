#pragma once

#include "ga.hpp"
#include "utility.hpp"

#include <fmt/color.h>
#include <fmt/format.h>
#include <fmt/ranges.h>

namespace fmt
{
template <typename T>
struct styled_object
{
    text_style ts;
    const T& value;
};

template <typename T>
[[nodiscard]] constexpr auto operator*(text_style lhs, const T& rhs) noexcept
{
    return styled_object<T>{lhs, rhs};
}

template <typename T>
struct formatter<styled_object<T>>
{
    template <typename PC>
    constexpr auto parse(PC& ctx)
    {
        return ctx.begin();
    }

    template <typename FC>
    constexpr auto format(const styled_object<T>& so, FC& ctx)
    {
        auto styled_text = ::fmt::vformat(so.ts, "{}", basic_format_args<FC>{format_arg_store<FC, T>{so.value}});
        return format_to(ctx.out(), "{}", styled_text);
    }
};

template <typename T, size_t... E>
struct formatter<gal::entity<T, E...>>
{
    using type = gal::entity<T, E...>;
    using elements = std::tuple<gal::element<E>...>;

    template <typename PC>
    constexpr auto parse(PC& ctx)
    {
        return ctx.begin();
    }

    template <typename FC>
    constexpr auto format(const type& entity, FC& ctx)
    {
        if constexpr (type::size == 0)
        {
            return ctx.out();
        }
        else
        {
            format_terms(entity, ctx, std::make_index_sequence<sizeof...(E)>{});
        }

        return ctx.out();
    }

    template <typename FC, size_t... Is>
    constexpr void format_terms(const type& entity, FC& ctx, std::index_sequence<Is...>)
    {
        (format_term(entity, std::integral_constant<size_t, Is>{}, ctx), ...);
    }

    template <size_t N, typename FC>
    constexpr void format_term(const type& entity, std::integral_constant<size_t, N>, FC& ctx)
    {
        if constexpr (N < sizeof...(E) - 1)
        {
            format_to(ctx.out(), "{}{} + ", entity.data[N], std::tuple_element_t<N, elements>{});
        }
        else
        {
            format_to(ctx.out(), "{} * {}", entity.data[N], std::tuple_element_t<N, elements>{});
        }
    }
};

template <int N, int D>
struct formatter<gal::rational<N, D>>
{
    using type = gal::rational<N, D>;

    template <typename PC>
    constexpr auto parse(PC& ctx)
    {
        return ctx.begin();
    }

    template <typename FC>
    constexpr auto format(const type&, FC& ctx)
    {
        if constexpr (N == 0)
        {
            format_to(ctx.out(), "0");
        }
        else if constexpr (D == 1)
        {
            format_to(ctx.out(), "{}", N);
        }
        else
        {
            format_to(ctx.out(), "{}/{}", N, D);
        }

        return ctx.out();
    }
};

template <typename Tag, typename Degree, size_t Order>
struct formatter<gal::generator<Tag, Degree, Order>>
{
    using type = gal::generator<Tag, Degree, Order>;

    template <typename PC>
    constexpr auto parse(PC& ctx)
    {
        return ctx.begin();
    }

    template <typename FC>
    constexpr auto format(const type&, FC& ctx)
    {
        // Untagged generators are displayed with a 0
        constexpr auto id = Tag::id == ~0ull ? 0 : Tag::id;
        if constexpr (Degree::value > 1)
        {
            return format_to(ctx.out(), "{}_{}^{}", fg(color::orange) * id, fg(color::chartreuse) * Tag::index, Degree::value);
        }
        else if constexpr (Degree::value == 1)
        {
            return format_to(ctx.out(), "{}_{}", fg(color::orange) * id, fg(color::chartreuse) * Tag::index);
        }
        else if constexpr (Degree::value == 0)
        {
            return format_to(ctx.out(), "{}_{}", fg(color::navajo_white) * 1);
        }
        else
        {
            return format_to(ctx.out(), "{}_{}^{{{}}}", fg(color::orange) * id, fg(color::chartreuse) * Tag::index, Degree::value);
        }
    }
};

template <typename Q, typename... Generators>
struct formatter<gal::monomial<Q, Generators...>>
{
    template <typename PC>
    constexpr auto parse(PC& ctx)
    {
        return ctx.begin();
    }

    template <typename FC>
    constexpr auto format(const gal::monomial<Q, Generators...>& addend, FC& ctx)
    {
        if constexpr (Q::is_zero)
        {
            return ctx.out();
        }
        else if constexpr (Q::num == -1 && Q::den == 1)
        {
            if constexpr (sizeof...(Generators) == 0)
            {
                format_to(ctx.out(), "-1");
            }
            else
            {
                format_to(ctx.out(), "-");
            }
        }
        else if constexpr (Q::num != 1 || Q::den != 1)
        {
            if constexpr (Q::den == 1)
            {
                const char* str = sizeof...(Generators) == 0 ? "{}" : "{} *";
                format_to(ctx.out(), str, fg(color::misty_rose) * Q::num);
            }
            else
            {
                const char* str = sizeof...(Generators) == 0 ? "{}/{}" : "{}/{} *";
                format_to(ctx.out(), str, fg(color::misty_rose) * Q::num, fg(color::misty_rose) * Q::den);
            }
        }
        else if constexpr (Q::num == 1 && Q::den == 1 && sizeof...(Generators) == 0)
        {
            format_to(ctx.out(), "{}", fg(color::misty_rose) * 1);
        }

        if constexpr (sizeof...(Generators) > 0)
        {
            format_generators(ctx, std::tuple<Generators...>{});
        }
        return ctx.out();
    }

    template <typename FC, typename G, typename... Gs>
    constexpr auto format_generators(FC& ctx, std::tuple<G, Gs...>)
    {
        if constexpr (sizeof...(Gs) > 0)
        {
            format_to(ctx.out(), "{}*", G{});
            format_generators(ctx, std::tuple<Gs...>{});
        }
        else
        {
            format_to(ctx.out(), "{}", G{});
        }
    }
};

template <size_t E>
struct formatter<gal::element<E>>
{
    template <typename PC>
    constexpr auto parse(PC& ctx)
    {
        return ctx.begin();
    }

    template <typename FC>
    constexpr auto format(const gal::element<E>& element, FC& ctx)
    {
        if constexpr (E == 0)
        {
            return ctx.out();
        }
        else
        {
            auto elem = E;
            while (gal::count_bits(elem) > 1)
            {
                auto index = gal::count_trailing_zeros(elem);
                format_to(ctx.out(), "e{}^", index);
                elem = (elem >> (index + 1)) << (index + 1);
            }

            return format_to(ctx.out(), "e{}", gal::count_trailing_zeros(elem));
        }
    }
};

template <typename E, typename... Monomials>
struct formatter<gal::term<E, Monomials...>>
{
    constexpr static text_style bstyle = fg(color::yellow);
    template <typename PC>
    constexpr auto parse(PC& ctx) const
    {
        return ctx.begin();
    }

    template <typename FC>
    constexpr auto format(const gal::term<E, Monomials...>& term, FC& ctx) const
    {
        if constexpr (sizeof...(Monomials) == 0)
        {
            return format_to(ctx.out(), "0");
        }

        format_to(ctx.out(), "[");
        if constexpr (sizeof...(Monomials) > 0)
        {
            format_addends(ctx, std::tuple<Monomials...>{});
        }
        if constexpr (E::value == 0)
        {
            format_to(ctx.out(), "]");
        }
        else
        {
            format_to(ctx.out(), "]{}", bstyle * E{});
        }
        return ctx.out();
    }

    template <typename FC, typename M, typename... Ms>
    constexpr auto format_addends(FC& ctx, std::tuple<M, Ms...>) const
    {
        if constexpr (sizeof...(Ms) > 0)
        {
            format_to(ctx.out(), "{} + ", M{});
            format_addends(ctx, std::tuple<Ms...>{});
        }
        else
        {
            format_to(ctx.out(), "{}", M{});
        }
    }
};

template <typename... Terms>
struct formatter<::gal::multivector<void, Terms...>>
{
    template <typename PC>
    constexpr auto parse(PC& ctx)
    {
        return ctx.begin();
    }

    template <typename FC>
    constexpr auto format(const gal::multivector<void, Terms...>& multivector, FC& ctx)
    {
        if constexpr (sizeof...(Terms) > 0)
        {
            format_terms(ctx, std::tuple<Terms...>{}, true);
        }
        else
        {
            format_to(ctx.out(), "0");
        }
        return ctx.out();
    }

    template <typename FC, typename T, typename... Ts>
    constexpr auto format_terms(FC& ctx, std::tuple<T, Ts...>, bool first)
    {
        if constexpr (sizeof...(Ts) > 0)
        {
            if constexpr (!T::is_zero)
            {
                if (!first)
                {
                    format_to(ctx.out(), " + {}", T{});
                }
                else
                {
                    format_to(ctx.out(), "{}", T{});
                }
                format_terms(ctx, std::tuple<Ts...>{}, false);
            }
            else
            {
                format_terms(ctx, std::tuple<Ts...>{}, first);
            }
        }
        else
        {
            if constexpr (!T::is_zero)
            {
                if (first)
                {
                    format_to(ctx.out(), "{}", T{});
                }
                else
                {
                    format_to(ctx.out(), " + {}", T{});
                }
            }
        }
    }
};
} // namespace fmt