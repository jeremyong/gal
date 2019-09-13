#pragma once

#include "finite_algebra.hpp"
#include "utility.hpp"

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <fmt/color.h>

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

    template <typename Degree, size_t ID>
    struct formatter<gal::factor<Degree, ID>>
    {
        template <typename PC>
        constexpr auto parse(PC& ctx)
        {
            return ctx.begin();
        }

        template <typename FC>
        constexpr auto format(const gal::factor<Degree, ID>& factor, FC& ctx)
        {
            if constexpr (Degree::value > 1)
            {
                return format_to(ctx.out(), "{}^{}", fg(color::orange) * ID, Degree::value);
            }
            else if constexpr (Degree::value == 1)
            {
                return format_to(ctx.out(), "{}", fg(color::orange) * ID);
            }
            else if constexpr (Degree::value == 0)
            {
                return format_to(ctx.out(), "{}", fg(color::navajo_white) * 1);
            }
            else
            {
                return format_to(ctx.out(), "{}^{{{}}}", fg(color::orange) * ID, Degree::value);
            }
        }
    };

    template <typename Multiplier, typename... Factors>
    struct formatter<gal::monomial<Multiplier, Factors...>>
    {
        template <typename PC>
        constexpr auto parse(PC& ctx)
        {
            return ctx.begin();
        }

        template <typename FC>
        constexpr auto format(const gal::monomial<Multiplier, Factors...>& addend, FC& ctx)
        {
            if constexpr (Multiplier::value == 0)
            {
                return ctx.out();
            }
            else if constexpr (Multiplier::value == -1)
            {
                format_to(ctx.out(), "-");
            }
            else if constexpr (Multiplier::value != 1)
            {
                format_to(ctx.out(), "{}*", fg(color::misty_rose) * Multiplier::value);
            }

            if constexpr (sizeof...(Factors) > 0)
            {
                format_factors(ctx, std::tuple<Factors...>{});
            }
            return ctx.out();
        }

        template <typename FC, typename F, typename... Fs>
        constexpr auto format_factors(FC& ctx, std::tuple<F, Fs...>)
        {
            if constexpr (sizeof...(Fs) > 0)
            {
                format_to(ctx.out(), "{}*", F{});
                format_factors(ctx, std::tuple<Fs...>{});
            }
            else
            {
                format_to(ctx.out(), "{}", F{});
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
        constexpr auto parse(PC& ctx)
        {
            return ctx.begin();
        }

        template <typename FC>
        constexpr auto format(const gal::term<E, Monomials...>& term, FC& ctx)
        {
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
        constexpr auto format_addends(FC& ctx, std::tuple<M, Ms...>)
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
    struct formatter<gal::multivector<void, Terms...>>
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
                format_terms(ctx, std::tuple<Terms...>{});
            }
            else
            {
                format_to(ctx.out(), "0");
            }
            return ctx.out();
        }

        template <typename FC, typename T, typename... Ts>
        constexpr auto format_terms(FC& ctx, std::tuple<T, Ts...>)
        {
            if constexpr (sizeof...(Ts) > 0)
            {
                format_to(ctx.out(), "{} + ", T{});
                format_terms(ctx, std::tuple<Ts...>{});
            }
            else
            {
                format_to(ctx.out(), "{}", T{});
            }
        }
    };
}