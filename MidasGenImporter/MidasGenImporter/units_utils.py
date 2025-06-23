
# TODO: check if all strings are compatible with midasgen (all uppercase, no spaces, etc.)

class unit_system:
    # Force conversions (base: newton)
    force_units = {
        # SI Units
        'N': 1.0,
        'KN': 1e3,
        'MN': 1e6,
        'dyne': 1e-5,
        'kgf': 9.80665,
        
        # US Customary
        'lbf': 4.44822,
        'kip': 4448.22,
        'poundal': 0.138255
    }

    # Length conversions (base: meter)
    length_units = {
        # SI Units
        'MM': 1e-3,
        'CM': 1e-2,
        'DM': 1e-1,
        'M': 1.0,
        'DAM': 10.0,
        'HM': 100.0,
        'KM': 1e3,

        # US Customary
        'IN': 0.0254,
        'FT': 0.3048,
        'YD': 0.9144,
        'MI': 1609.344,
        'MIL': 2.54e-5,
        'THOU': 2.54e-5
    }

    @staticmethod
    def F(value: float, from_unit: str, to_unit: str) -> float:
        if from_unit not in unit_system.force_units or to_unit not in unit_system.force_units:
            raise ValueError(f"Unsupported force unit: {from_unit} or {to_unit}")
        return value * unit_system.force_units[from_unit] / unit_system.force_units[to_unit]

    @staticmethod
    def L(value: float, from_unit: str, to_unit: str) -> float:
        if from_unit not in unit_system.length_units or to_unit not in unit_system.length_units:
            raise ValueError(f"Unsupported length unit: {from_unit} or {to_unit}")
        return value * unit_system.length_units[from_unit] / unit_system.length_units[to_unit]

    @staticmethod
    def T(value: float, from_unit: str, to_unit: str) -> float:
        from_unit = from_unit.upper()
        to_unit = to_unit.upper()

        if from_unit == to_unit:
            return value

        # Convert to Celsius
        if from_unit == 'C':
            c = value
        elif from_unit == 'F':
            c = (value - 32) * 5 / 9
        elif from_unit == 'K':
            c = value - 273.15
        elif from_unit == 'R':
            c = (value - 491.67) * 5 / 9
        else:
            raise ValueError(f"Unsupported temperature unit: {from_unit}")

        # Convert from Celsius to target unit
        if to_unit == 'C':
            return c
        elif to_unit == 'F':
            return c * 9 / 5 + 32
        elif to_unit == 'K':
            return c + 273.15
        elif to_unit == 'R':
            return (c + 273.15) * 9 / 5
        else:
            raise ValueError(f"Unsupported temperature unit: {to_unit}")