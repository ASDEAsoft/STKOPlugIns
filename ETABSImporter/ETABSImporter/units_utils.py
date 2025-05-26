
class unit_system:
    # Force conversions (base: newton)
    force_units = {
        # SI Units
        'N': 1.0,
        'kN': 1e3,
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
        'mm': 1e-3,
        'cm': 1e-2,
        'dm': 1e-1,
        'm': 1.0,
        'dam': 10.0,
        'hm': 100.0,
        'km': 1e3,

        # US Customary
        'in': 0.0254,
        'ft': 0.3048,
        'yd': 0.9144,
        'mi': 1609.344,
        'mil': 2.54e-5,
        'thou': 2.54e-5
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