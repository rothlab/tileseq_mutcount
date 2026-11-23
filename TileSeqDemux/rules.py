"""Small objects that say what can be wrong with one row of a table, so that a new check is one new Rule.

    Rule      one thing that can be wrong: its name, the message to show, and the test that finds it
    RuleSet   the ordered rules for one kind of row, and the code that applies them
    Verdict   what applying a RuleSet to one row found

A row class gives a RuleSet what it needs by providing `details(context)`, a dict of the values the messages can use.
The messages are written with {placeholders}, and a RuleSet is told which placeholders the row provides, so a typo in
a message or a repeated rule name is an error when the module loads, not when a rare bad row finally turns up.
Standard library only.
"""
import re
import string
from dataclasses import dataclass
from types import SimpleNamespace
from typing import Callable, List, NamedTuple

BLOCKS = "blocks"
REPORTS = "reports"
WARNS = "warns"
EFFECTS = (BLOCKS, REPORTS, WARNS)


@dataclass(frozen=True)
class Rule:
    """One thing that can be wrong with a row.

    `broken_when(row, context)` returns True when the row breaks the rule."""
    name: str
    message: str
    broken_when: Callable
    effect: str = BLOCKS


class Verdict(NamedTuple):
    names: List[str]          # the rules the row broke, in order
    problems: List[str]       # their messages when the effect is BLOCKS or REPORTS
    warnings: List[str]       # their messages when the effect is WARNS
    blocked: bool             # True when the row must be left out or cannot be used

class Facts(SimpleNamespace):
    """The values a RuleSet needs, as keyword arguments: Facts(count=2, label='CDS end:')."""

    def details(self, context=None):
        return dict(vars(self))

class RuleSet:
    """The ordered rules for one kind of thing. Rules are tried in order and the effect of a broken rule decides what follows."""

    def __init__(self, rules, details):
        self.rules = tuple(rules)
        self._details = frozenset(details)
        self._validate()

    def check(self, subject, context=None):
        """Apply the rules to `subject`. `context` is whatever the rules need to know about the other rows."""
        details = subject.details(context)
        missing = self._details - set(details)
        if missing:
            raise KeyError(f"the thing being checked does not provide {sorted(missing)}")
        names, problems,warnings = [], [], []
        for rule in self.rules:
            if not rule.broken_when(subject, context):
                continue
            names.append(rule.name)
            (warnings if rule.effect == WARNS else problems).append(rule.message.format(**details))
            if rule.effect == BLOCKS:
                return Verdict(names, problems, warnings, blocked=True)
        return Verdict(names, problems, warnings, blocked=False)

    def _validate(self):
        names = [rule.name for rule in self.rules]
        repeated = sorted({name for name in names if names.count(name) > 1})
        if repeated:
            raise ValueError(f"rule names must be unique, repeated: {repeated}")
        for rule in self.rules:
            if rule.effect not in EFFECTS:
                raise ValueError(f"rule {rule.name!r} has effect {rule.effect!r}; it must be one of {EFFECTS}")
            unknown = self._placeholders(rule.message) - self._details
            if unknown:
                raise ValueError(f"rule {rule.name!r} uses {sorted(unknown)}, which the thing being checked does not provide")

    @staticmethod
    def _placeholders(message):
        fields = (field for _text, field, _spec, _conversion in string.Formatter().parse(message) if field)
        return {re.split(r"[.\[]", field)[0] for field in fields}