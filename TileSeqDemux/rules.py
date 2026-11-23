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
from typing import Callable, List, NamedTuple


@dataclass(frozen=True)
class Rule:
    """One thing that can be wrong with a row.

    `broken_when(row, context)` returns True when the row breaks the rule. `blocks_row=False` means the problem is
    reported but the row is kept."""
    name: str
    message: str
    broken_when: Callable
    blocks_row: bool = True


class Verdict(NamedTuple):
    names: List[str]          # the rules the row broke, in order
    messages: List[str]       # their messages, filled in
    blocked: bool             # True when the row must be left out


class RuleSet:
    """The ordered rules for one kind of row.

    Rules are tried in order. A broken rule that blocks the row ends the checking, so a later rule can rely on the
    earlier ones (a range rule can assume the numbers were read). A broken rule that does not block records its
    message and the checking goes on."""

    def __init__(self, rules, details):
        self.rules = tuple(rules)
        self._details = frozenset(details)
        self._validate()

    def check(self, row, context):
        """Apply the rules to `row`. `context` is whatever the rules need to know about the other rows."""
        details = row.details(context)
        names, messages = [], []
        for rule in self.rules:
            if rule.broken_when(row, context):
                names.append(rule.name)
                messages.append(rule.message.format(**details))
                if rule.blocks_row:
                    return Verdict(names, messages, True)
        return Verdict(names, messages, False)

    def _validate(self):
        names = [rule.name for rule in self.rules]
        repeated = sorted({name for name in names if names.count(name) > 1})
        if repeated:
            raise ValueError(f"rule names must be unique, repeated: {repeated}")
        for rule in self.rules:
            unknown = self._placeholders(rule.message) - self._details
            if unknown:
                raise ValueError(f"rule {rule.name!r} uses {sorted(unknown)}, which the row does not provide")

    @staticmethod
    def _placeholders(message):
        fields = (field for _text, field, _spec, _conversion in string.Formatter().parse(message) if field)
        return {re.split(r"[.\[]", field)[0] for field in fields}