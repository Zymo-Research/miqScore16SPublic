defaultAllowableUniqueIDCharacters = "."

def isAlphaNumeric(character:str, allowedCharacters:str = defaultAllowableUniqueIDCharacters, replacement = None):
    if not type(character) == str:
        return False
    if not character:
        return False
    if len(character) > 1:
        return alphaNumericString(character, allowedCharacters, replacement)
    if allowedCharacters:
        if character in allowedCharacters:
            return True
    if character.isalnum():
        return True
    return False


def alphaNumericString(string:str, allowedCharacters:str = defaultAllowableUniqueIDCharacters, replacement = None):
    if not string:
        return False
    validatedString = ""
    for character in string:
        if not isAlphaNumeric(character, allowedCharacters):
            if replacement != None:
                validatedString += replacement
            else:
                return False
        else:
            validatedString += character
    return validatedString
