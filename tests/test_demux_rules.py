# These tests were generated with Claude (Anthropic) and then reviewed, edited and run by the maintainers.
"""
Tests for the Rule, RuleSet, Verdict and Facts classes (ticket 1.3): rules are tried in order, a rule's effect decides
whether the checking stops, a problem is kept apart from a warning, and mistakes in a rule table are errors as soon as
they are seen.

Run with:  python -m pytest tests/test_demux_rules.py -v     (from the repository root)
Without pytest: python -m unittest discover -s tests -p "test_demux_*.py"
"""
import dataclasses
import unittest

from TileSeqDemux.rules import BLOCKS, REPORTS, WARNS, Facts, Rule, RuleSet, Verdict


class Row:
    def __init__(self, value):
        self.value = value

    def details(self, context):
        return {"value": self.value, "limit": context}


def make_rules(*rules):
    return RuleSet(rules, details=("value", "limit"))


TOO_BIG = Rule("too_big", "{value} is above {limit}", lambda row, limit: row.value > limit)
NEGATIVE = Rule("negative", "{value} is negative", lambda row, limit: row.value < 0)
ODD = Rule("odd", "{value} is odd", lambda row, limit: row.value % 2 == 1, effect=REPORTS)
LONG = Rule("long", "{value} has more than one digit", lambda row, limit: row.value >= 10, effect=WARNS)


class CheckTests(unittest.TestCase):
    def test_a_row_that_breaks_nothing(self):
        self.assertEqual(make_rules(TOO_BIG, NEGATIVE).check(Row(4), 10), Verdict([], [], [], False))

    def test_a_broken_rule_gives_its_name_and_message(self):
        verdict = make_rules(TOO_BIG, NEGATIVE).check(Row(40), 10)
        self.assertEqual(verdict.names, ["too_big"])
        self.assertEqual(verdict.problems, ["40 is above 10"])
        self.assertEqual((verdict.warnings, verdict.blocked), ([], True))

    def test_a_blocking_rule_stops_the_checking(self):
        """Both rules would fire on -3 with a limit of -5, but only the first one is reported."""
        self.assertEqual(make_rules(TOO_BIG, NEGATIVE).check(Row(-3), -5).names, ["too_big"])

    def test_order_decides_which_rule_is_reported(self):
        self.assertEqual(make_rules(NEGATIVE, TOO_BIG).check(Row(-3), -5).names, ["negative"])

    def test_a_later_rule_can_rely_on_an_earlier_blocking_rule(self):
        """The second rule would raise TypeError on None, so it must never be called."""
        parsed = Rule("not_a_number", "{value!r} is not a number", lambda row, limit: row.value is None)
        compared = Rule("too_big", "{value} is above {limit}", lambda row, limit: row.value > limit)
        self.assertEqual(make_rules(parsed, compared).check(Row(None), 10).names, ["not_a_number"])

    def test_a_reporting_rule_records_a_problem_and_the_checking_goes_on(self):
        verdict = make_rules(ODD, TOO_BIG).check(Row(11), 10)
        self.assertEqual(verdict.names, ["odd", "too_big"])
        self.assertEqual(verdict.problems, ["11 is odd", "11 is above 10"])
        self.assertTrue(verdict.blocked)

    def test_a_reporting_rule_alone_does_not_block(self):
        verdict = make_rules(TOO_BIG, ODD).check(Row(7), 10)
        self.assertEqual((verdict.names, verdict.problems, verdict.blocked), (["odd"], ["7 is odd"], False))

    def test_a_warning_is_kept_apart_from_the_problems_and_never_blocks(self):
        verdict = make_rules(LONG, TOO_BIG).check(Row(12), 20)
        self.assertEqual((verdict.names, verdict.problems, verdict.warnings, verdict.blocked),
                         (["long"], [], ["12 has more than one digit"], False))

    def test_a_warning_and_a_problem_can_come_from_one_row(self):
        verdict = make_rules(LONG, TOO_BIG).check(Row(12), 5)
        self.assertEqual((verdict.problems, verdict.warnings, verdict.blocked), (["12 is above 5"], ["12 has more than one digit"], True))

    def test_the_context_reaches_the_rules_and_the_messages(self):
        self.assertEqual(make_rules(TOO_BIG).check(Row(5), 3).problems, ["5 is above 3"])

    def test_a_rule_cannot_be_changed_after_it_is_made(self):
        with self.assertRaises(dataclasses.FrozenInstanceError):
            TOO_BIG.name = "other"

    def test_a_rule_blocks_unless_told_otherwise(self):
        self.assertEqual(TOO_BIG.effect, BLOCKS)


class FactsTests(unittest.TestCase):
    def test_values_are_read_as_attributes_and_given_to_the_messages(self):
        rules = RuleSet([Rule("empty", "{what} has {count} rows", lambda facts, context: facts.count == 0)], details=("what", "count"))
        verdict = rules.check(Facts(what="the table", count=0))
        self.assertEqual(verdict.problems, ["the table has 0 rows"])

    def test_the_context_is_optional(self):
        rules = RuleSet([Rule("never", "x", lambda facts, context: False)], details=())
        self.assertEqual(rules.check(Facts(a=1)), Verdict([], [], [], False))

    def test_a_value_the_rule_set_expects_but_was_not_given_is_an_error_even_if_no_rule_fires(self):
        rules = RuleSet([Rule("never", "{count}", lambda facts, context: False)], details=("count",))
        with self.assertRaisesRegex(KeyError, "count"):
            rules.check(Facts(other=1))


class MistakesAreCaughtWhenTheTableIsBuilt(unittest.TestCase):
    def test_a_repeated_rule_name(self):
        with self.assertRaisesRegex(ValueError, "unique.*too_big"):
            make_rules(TOO_BIG, Rule("too_big", "again", lambda row, limit: False))

    def test_a_placeholder_the_row_does_not_provide(self):
        with self.assertRaisesRegex(ValueError, r"'typo'.*\['valu'\]"):
            make_rules(Rule("typo", "{valu} is wrong", lambda row, limit: True))

    def test_an_effect_that_does_not_exist(self):
        with self.assertRaisesRegex(ValueError, "'odd'.*'sometimes'"):
            make_rules(Rule("odd", "{value}", lambda row, limit: True, effect="sometimes"))

    def test_conversions_and_format_specs_are_understood(self):
        make_rules(Rule("fine", "{value!r} {limit:>5}", lambda row, limit: False))

    def test_attribute_and_index_placeholders_are_checked_by_their_first_part(self):
        make_rules(Rule("fine", "{value.real} {limit[0]}", lambda row, limit: False))
        with self.assertRaises(ValueError):
            make_rules(Rule("typo", "{other.real}", lambda row, limit: False))

    def test_doubled_braces_are_not_placeholders(self):
        make_rules(Rule("fine", "a literal {{brace}} and {value}", lambda row, limit: False))


if __name__ == "__main__":
    unittest.main()