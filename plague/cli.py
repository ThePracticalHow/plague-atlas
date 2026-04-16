"""
plague CLI — pathogen detection from the command line.

    plague detect data.h5ad
    plague detect data.h5ad --group cell_type -o results.json
    plague food garlic
    plague food --report
    plague panels
    plague version
"""

import argparse
import sys


def cmd_detect(args):
    from .detect import detect
    out = args.output or args.path.replace('.h5ad', '_plague_detect.json')
    detect(args.path, group_key=args.group, output_path=out)


def cmd_food(args):
    from .food import FoodAtlas
    atlas = FoodAtlas()

    if args.report:
        atlas.war_report()
        return

    if not args.name:
        atlas.war_report()
        return

    food = ' '.join(args.name)
    result = atlas.score_food(food)
    print(f'\n  Food: {food}')
    print(f'  Score: {result["score"]:+d}')
    if result.get('unknown'):
        print(f'  (not in database)')
    else:
        if result['antifungal']:
            print(f'  Antifungal compounds:')
            for name, potency in result['antifungal']:
                print(f'    {name}: potency {potency}')
        if result['profungal']:
            print(f'  Profungal compounds:')
            for name, potency in result['profungal']:
                print(f'    {name}: potency {potency}')


def cmd_panels(args):
    from .detect import ALARM_GENES, SPORE_GENES, MELANIN_GENES

    print('\nPlague Detection Panels')
    print('=' * 60)

    print(f'\n  ALARM panel ({len(ALARM_GENES)} genes) — fungal detection:')
    for gene, desc in ALARM_GENES.items():
        print(f'    {gene:10s}  {desc}')

    print(f'\n  SPORE panel ({len(SPORE_GENES)} genes) — lysosomal expansion:')
    for gene, desc in list(SPORE_GENES.items())[:10]:
        print(f'    {gene:10s}  {desc}')
    if len(SPORE_GENES) > 10:
        print(f'    ... and {len(SPORE_GENES) - 10} more')

    print(f'\n  MELANIN panel ({len(MELANIN_GENES)} genes) — melanin machinery:')
    for gene, desc in MELANIN_GENES.items():
        print(f'    {gene:10s}  {desc}')


def cmd_tensor_detect(args):
    """Combined: compute coupling tensor + plague detection."""
    from coherence.tensor import run as coherence_run
    from .detect import detect

    print("Phase 1: Coupling tensor")
    tensor_out = args.path.replace('.h5ad', '_tensor.json')
    coherence_run(args.path, output_path=tensor_out, group_key=args.group)

    print("\nPhase 2: Plague detection")
    detect_out = args.path.replace('.h5ad', '_plague_detect.json')
    detect(args.path, group_key=args.group, output_path=detect_out)

    print(f"\nTensor:  {tensor_out}")
    print(f"Detect:  {detect_out}")


def cmd_version(args):
    from . import __version__
    print(f'plague {__version__}')


def main():
    parser = argparse.ArgumentParser(
        prog='plague',
        description='Pathogen detection engine. Fungal scanner, diagnosis, dietary war scoring.',
    )
    sub = parser.add_subparsers(dest='command')

    p_detect = sub.add_parser('detect', help='Detect pathogen indicators in expression data')
    p_detect.add_argument('path', help='Path to .h5ad file')
    p_detect.add_argument('--group', '-g', help='Column in obs to group by')
    p_detect.add_argument('--output', '-o', help='Output JSON path')
    p_detect.set_defaults(func=cmd_detect)

    p_food = sub.add_parser('food', help='Score a food or print the war report')
    p_food.add_argument('name', nargs='*', help='Food name to score')
    p_food.add_argument('--report', action='store_true', help='Print full war report')
    p_food.set_defaults(func=cmd_food)

    p_panels = sub.add_parser('panels', help='Show detection gene panels')
    p_panels.set_defaults(func=cmd_panels)

    p_full = sub.add_parser('full', help='Combined: coupling tensor + plague detection')
    p_full.add_argument('path', help='Path to .h5ad file')
    p_full.add_argument('--group', '-g', help='Group column')
    p_full.set_defaults(func=cmd_tensor_detect)

    p_ver = sub.add_parser('version', help='Print version')
    p_ver.set_defaults(func=cmd_version)

    args = parser.parse_args()
    if not args.command:
        parser.print_help()
        return

    args.func(args)


if __name__ == '__main__':
    main()
