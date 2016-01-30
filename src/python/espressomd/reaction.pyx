include "myconfig.pxi"
cimport reaction
cimport globals
cimport utils

_reaction_is_instanciated = False

IF CATALYTIC_REACTIONS:
    cdef class reaction:
        def __init__(self, product_type=None, reactant_type=None, catalyzer_type=None, ct_range=None,\
                     ct_rate=None, eq_rate=0.0, react_once=False, swap=False):
            """
            Initialize the reaction.  Keep in mind, that there may be
            only one reaction enabled.  You will be warned if you
            still instanciate a second one.

            Dictionary initialization is also feasible.  With a
            dictionary named `params' it reads

                espressomd.reaction.reaction(params)
            """

            global _reaction_is_instanciated
            if _reaction_is_instanciated:
                print "WARNING: There can be only one reaction! I'm overwriting any previous ones."
            _reaction_is_instanciated = True

            if isinstance(product_type, dict):
                params = product_type.copy()
                product_type = params.get('product_type')
                reactant_type = params.get('reactant_type')
                catalyzer_type = params.get('catalyzer_type')
                ct_range = params.get('ct_range')
                ct_rate = params.get('ct_rate')
                eq_rate = params.get('eq_rate')
                react_once = bool(params.get('react_once'))
                swap = bool(params.get('swap'))

            utils.check_type_or_throw_except(product_type,   1, int,   "product_type has to be an int")
            utils.check_type_or_throw_except(reactant_type,  1, int,   "reactant_type has to be an int")
            utils.check_type_or_throw_except(catalyzer_type, 1, int,   "catalyzer_type has to be an int")
            utils.check_type_or_throw_except(ct_range,       1, float, "ct_range has to be a float")
            utils.check_type_or_throw_except(ct_rate,        1, float, "ct_rate has to be a float")
            utils.check_type_or_throw_except(eq_rate,        1, float, "eq_rate has to be a float")
            utils.check_type_or_throw_except(react_once,     1, bool,  "react_once has to be a bool")
            utils.check_type_or_throw_except(swap,           1, bool,  "swap has to be a bool")

            if( product_type == reactant_type or product_type == catalyzer_type or catalyzer_type == reactant_type ):
                raise Exception("One particle type cannot be a part more than one reaction species!")
            globals.reaction.product_type   = product_type
            globals.reaction.reactant_type  = reactant_type
            globals.reaction.catalyzer_type = catalyzer_type

            if (ct_range < 0):
                raise ValueError("Negative equilibrium reaction rate contstant is not allowed!")
            globals.reaction.range = ct_range

            if (ct_rate < 0):
                raise ValueError("Negative catalytic reaction rate constant is not allowed!")
            globals.reaction.ct_rate = ct_rate

            if (eq_rate < 0 and abs(reaction.eq_rate + 1.0) > 0.001 ):
                raise ValueError("Negative equilibrium reaction rate contstant is not allowed!")
            globals.reaction.eq_rate = eq_rate

            globals.reaction.sing_mult = int(react_once)
            globals.reaction.swap = int(swap)

            mpi_setup_reaction()


        def inhibit(self):
            """
            Inhibit the reaction, i.e. set the reaction rate to 0.0
            """
            globals.reaction.ct_rate = 0.0
            mpi_setup_reaction()

        def get_params(self):
            """
            Returns a dictionary with the reaction parameters.  The
            format is suitable for dictionary initialization of the
            reaction.
            """
            return {
                'reactant_type': globals.reaction.reactant_type,
                'catalyzer_type': globals.reaction.catalyzer_type,
                'product_type': globals.reaction.product_type,
                'ct_range': globals.reaction.range,
                'ct_rate': globals.reaction.ct_rate,
                'eq_rate': globals.reaction.eq_rate,
                'react_once': bool(globals.reaction.sing_mult),
                'swap': bool(globals.reaction.swap)
            }

        def __repr__(self):
            return "espressomd.reaction()"


        def __str__(self):
            return str(self.get_params())
