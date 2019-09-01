#Numerical Validations

def isInteger(value):
    return type(value) == int


def isNumber(value):
    return type(value) in [int, float]


def isPositive(value):
    if isNumber(value):
        return value > 0
    else:
        return False


def isNotNegative(value):
    if isNumber(value):
        return value >= 0
    else:
        return False


def isPositiveInteger(value):
    return type(value) == int and isPositive(value)


def isNonNegativeInteger(value):
    return type(value) == int and isNotNegative(value)

