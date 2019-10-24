-- Copyright (C) 2017 Henri Menke
--
-- Distributed under the Boost Software License, Version 1.0. (See
-- accompanying file LICENSE_1_0.txt or copy at
-- http://www.boost.org/LICENSE_1_0.txt)

local lpeg = require"lpeg"
local C, P, R, S, V = lpeg.C, lpeg.P, lpeg.R, lpeg.S, lpeg.V

local white = S(" \t") ^ 0

local digit    = R("09")
local dot      = P(".")
local eE       = S("eE")
local sign     = S("+-")^-1
local mantissa = digit^1 * dot * digit^0 + dot * digit^1 + digit^1
local exponent = (eE * sign * digit^1)^-1
local real     = white * mantissa * exponent * white / tonumber
local power    = white * C(P("^")) * white
local muldiv   = white * C(S("/*%")) * white
local addsub   = white * C(S("+-")) * white
local open     = white * P("(") * white
local close    = white * P(")") * white
local comma    = white * P(",") * white
local variable = white * C(R("az","AZ") * (R("az","AZ","09") + P("_"))^0) * white

-- Lookup for functions
local ufunc = {
    ["+"] = function(x) return x end, ["-"] = function(x) return -x end,
    sgn = function(x) return (0<x and 1 or 0) - (x<0 and 1 or 0) end,
    isinf = function(x) return (x == math.huge or x == -math.huge) and 1 or 0 end,
    isnan = function(x) return x ~= x and 1 or 0 end,
    abs  = math.abs,  acos  = math.acos,  asin = math.asin, atan  = math.atan,
    ceil = math.ceil, cos   = math.cos,   cosh = math.cosh, deg   = math.deg,
    exp  = math.exp,  floor = math.floor, log  = math.log,  log10 = math.log10,
    rad  = math.rad,  sin   = math.sin,   sinh = math.sinh, sqrt  = math.sqrt,
    tan  = math.tan,  tanh  = math.tanh
}
local bfunc = {
    ["+"] = function(a,b) return a + b end,
    ["-"] = function(a,b) return a - b end,
    ["*"] = function(a,b) return a * b end,
    ["/"] = function(a,b) return a / b end,
    ["%"] = function(a,b) return a % b end,
    ["^"] = function(a,b) return a ^ b end,
    atan2 = math.atan2, max = math.max, min = math.min, pow = math.pow
}

-- Lookup for contants
local const = {
    e = math.exp(1), inf = math.huge, nan = 0/0, phi = (1+math.sqrt(5))/2, pi = math.pi
}

-- Generate rules from lookup tables
local function generate(lookup)
    local rule
    for k,_ in pairs(lookup) do
        rule = rule and rule + C(P(k)) or C(P(k))
    end
    return rule
end

local ufunc_rule = generate(ufunc)
local bfunc_rule = generate(bfunc)
local const_rule = generate(const)

-- Evaluate AST recursively
local function eval(t,s)
    if type(t) == "table" then
        if t.tag == "u" then
            return ufunc[t.op](eval(t.right,s))
        elseif t.tag == "b" then
            return bfunc[t.op](eval(t.left,s),eval(t.right,s))
        else
            error("Cannot happen")
        end
    elseif type(t) == "number" then
        return t
    elseif const[t] then
        return const[t]
    elseif s and s[t] then
        return s[t]
    else
        error(string.format("Unknown operand '%s'",t))
    end
end

-- Insert binary node into AST
local function binary(rule)
    local function recurse(left,op,right,...)
        if op then
            return recurse({ tag = "b", op = op, left = left, right = right },...)
        else
            return left
        end
    end
    return rule / recurse
end

-- Insert unary node into AST
local function unary(rule)
    return rule / function(op,right)
        return { tag = "u", op = op, right = right }
    end
end

local grammar = P({
        "input",
        input      = V("expression") * -1,
        expression = binary( V("term") * ( addsub * V("term") )^0 ),
        term       = binary( V("factor") * ( muldiv * V("factor"))^0 ),
        factor     = binary( V("primary") * ( power * V("factor"))^0 ),
        primary    = real
            + open * V("expression") * close
            + unary( addsub * V("primary") )
            + unary( ufunc_rule * open * V("expression") * close )
            + binary( bfunc_rule * open * V("expression") * comma * V("expression") * close /
                      function(op,left,right) return left,op,right end ) -- switch around operands
            + const_rule
            + variable
})

local matheval = {
    parse = function(expr,symbols)
        return eval(assert(grammar:match(expr)),symbols)
    end
}

return matheval

--[[
print(matheval.parse("nan"), 0/0)
print(matheval.parse("inf"), math.huge)
print(matheval.parse("isinf(1e1000)"), ufunc.isinf(1e1000))
print(matheval.parse("isnan(0/0)"), ufunc.isnan(0/0))
print(matheval.parse("sgn(-5)"), ufunc.sgn(-5))
print(matheval.parse("6%4"), 6%4)
print(matheval.parse("3*(5-7^8)"), 3*(5-7^8))
print(matheval.parse("5 + 2*3 - 1 + 7 * 8"), 5 + 2*3 - 1 + 7 * 8)
print(matheval.parse("67 + 2 * 3 - 67 + 2/1 - 7"), 67 + 2 * 3 - 67 + 2/1 - 7)
print(matheval.parse("(5*7/5) + (23) - 5 * (98-4)/(6*7-42)"), (5*7/5) + (23) - 5 * (98-4)/(6*7-42))
print(matheval.parse("(2) + (17*2-30) * (5)+2 - (8/2)*4"), (2) + (17*2-30) * (5)+2 - (8/2)*4)
print(matheval.parse("2*3*4/8 -   5/2*4 +  6 + 0/3"), 2*3*4/8 -   5/2*4 +  6 + 0/3)
print(matheval.parse("2*3 - 4*5 + 6/3"), 2*3 - 4*5 + 6/3)
print(matheval.parse("tan(.7438E+0)"), math.tan(.7438E+0))
print(matheval.parse("sin(pi)"), math.sin(math.pi))
print(matheval.parse("pow(2,3)"), 2^3)
print(matheval.parse("min(2,3)"), math.min(2,3))
print(matheval.parse("max(2,3)"), math.max(2,3))
print(matheval.parse("y0_abc_32_def + 2",{ y0_abc_32_def = 1 }), 3)
--]]
