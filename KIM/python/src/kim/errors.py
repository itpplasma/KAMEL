"""Public exception hierarchy for the KIM Python interface."""


class KimError(Exception):
    """Base class for errors raised by the KIM Python interface."""


class ConfigurationError(KimError):
    """Raised when a KIM configuration cannot be constructed or loaded."""


class ProfileError(KimError):
    """Raised when radial profile input is missing, malformed, or inconsistent."""


class ExecutableError(KimError):
    """Raised when the KIM scientific executable cannot be resolved or inspected."""


class RunError(KimError):
    """Raised when a run repository or manifest operation fails."""
