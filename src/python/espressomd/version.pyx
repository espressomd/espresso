def major():
    """Prints the major version of Espresso.
    """
    return ESPRESSO_VERSION // 100000


def minor():
    """Prints the minor version of Espresso.
    """
    return (ESPRESSO_VERSION // 100) % 1000


def friendly():
    """Dot version of the version.
    """
    return "{}.{}".format(major(), minor())


def git_branch():
    """Git branch of the build if known, otherwise
       empty.
    """
    return GIT_BRANCH


def git_commit():
    """Git commit of the build if known, otherwise
       empty.
    """
    return GIT_COMMIT_HASH


def git_state():
    """Git state of the build if known, otherwise
       empty. State is "CLEAN" if the repository
       was not changed from git_commit(), "DIRTY"
       otherwise.
    """
    return GIT_STATE
