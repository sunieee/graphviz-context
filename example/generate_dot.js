let dot = `
digraph G {
    year1988 [label="1988"]
    year1989 [label="1989"]
    year1990 [label="1990"]
    year1991 [label="1991"]
    year1992 [label="1992"]
    year1993 [label="1993"]
    year1994 [label="1994"]
    year1995 [label="1995"]
    year1996 [label="1996"]
    year1997 [label="1997"]
    year1998 [label="1998"]
    year1999 [label="1999"]
    year2000 [label="2000"]
    year2001 [label="2001"]
    year2002 [label="2002"]
    year2003 [label="2003"]
    year2004 [label="2004"]
    year2005 [label="2005"]
    year2006 [label="2006"]
    year2007 [label="2007"]
    year2008 [label="2008"]
    year2009 [label="2009"]
    year2010 [label="2010"]
    year2011 [label="2011"]
    year2012 [label="2012"]
    year2013 [label="2013"]
    year2014 [label="2014"]
    year2015 [label="2015"]
    year2016 [label="2016"]
    year2017 [label="2017"]
    year2018 [label="2018"]
    year2019 [label="2019"]
    year2020 [label="2020"]
    year2021 [label="2021"]
    year2022 [label="2022"]
    1710476689 [label="2352"]
    1786904711 [label="64"]
    1843615162 [label="17"]
    2025768430 [label="4424"]
    2099471712 [label="22672"]
    2123284177 [label="29"]
    2134842679 [label="167"]
    2159528849 [label="15"]
    2421157170 [label="6"]
    2585630030 [label="234"]
    2592298275 [label="100"]
    2593383075 [label="141"]
    2751842161 [label="140"]
    2770487696 [label="5"]
    2950320139 [label="58"]
    2953243993 [label="1"]
    2953267151 [label="163"]
    2953501176 [label="4"]
    2963622136 [label="11"]
    2963865839 [label="98"]
    2970006822 [label="249"]
    2971040589 [label="90"]
    2980709326 [label="12"]
    3001314098 [label="0"]
    3011483492 [label="13"]
    3012026596 [label="5"]
    3093990297 [label="4"]
    { rank=same year1988 }
    { rank=same year1989 }
    { rank=same year1990 }
    { rank=same year1991 }
    { rank=same year1992 }
    { rank=same year1993 }
    { rank=same year1994 }
    { rank=same year1995 }
    { rank=same year1996 }
    { rank=same year1997 }
    { rank=same year1998 }
    { rank=same year1999 }
    { rank=same year2000 }
    { rank=same year2001 }
    { rank=same year2002 }
    { rank=same year2003 }
    { rank=same year2004 }
    { rank=same year2005 }
    { rank=same year2006 }
    { rank=same year2007 }
    { rank=same year2008 2025768430 }
    { rank=same year2009 }
    { rank=same year2010 }
    { rank=same year2011 }
    { rank=same year2012 1786904711 2123284177 2950320139 }
    { rank=same year2013 1843615162 2134842679 2953267151 }
    { rank=same year2014 1710476689 2099471712 }
    { rank=same year2015 2421157170 }
    { rank=same year2016 2159528849 2585630030 2963622136 2963865839 }
    { rank=same year2017 2592298275 2593383075 2751842161 2770487696 }
    { rank=same year2018 }
    { rank=same year2019 2953243993 2953501176 2970006822 2971040589 2980709326 }
    { rank=same year2020 3001314098 3011483492 3012026596 3093990297 }
    { rank=same year2021 }
    { rank=same year2022 }
    year1988->year1989
    year1989->year1990
    year1990->year1991
    year1991->year1992
    year1992->year1993
    year1993->year1994
    year1994->year1995
    year1995->year1996
    year1996->year1997
    year1997->year1998
    year1998->year1999
    year1999->year2000
    year2000->year2001
    year2001->year2002
    year2002->year2003
    year2003->year2004
    year2004->year2005
    year2005->year2006
    year2006->year2007
    year2007->year2008
    year2008->year2009
    year2009->year2010
    year2010->year2011
    year2011->year2012
    year2012->year2013
    year2013->year2014
    year2014->year2015
    year2015->year2016
    year2016->year2017
    year2017->year2018
    year2018->year2019
    year2019->year2020
    year2020->year2021
    year2021->year2022
    2025768430->2099471712
    2025768430->2159528849
    2025768430->2953501176
    2123284177->2134842679
    2950320139->2099471712
    2950320139->2134842679
    2950320139->2159528849
    2950320139->2953267151
    2134842679->2099471712
    2134842679->2159528849
    2134842679->2953501176
    2953267151->1843615162
    1843615162->2134842679
    1843615162->2953267151
    2099471712->2585630030
    2099471712->2592298275
    2099471712->2593383075
    2099471712->2980709326
    2099471712->3011483492
    2099471712->3012026596
    2592298275->2593383075
    2970006822->3093990297
}
`

let context_edges = {
    "l2010->2421157170": 1,
    "l2010->2963622136": 1,
    "l2006->2159528849": 1,
    "l2012->2134842679": 1,
    "l2012->2953267151": 1,
    "l2011->2123284177": 2,
    "l2011->2134842679": 1,
    "l2011->2953267151": 1,
    "l2011->1843615162": 1,
    "l2011->2099471712": 2,
    "l2011->2159528849": 1,
    "l2010->1843615162": 1,
    "2025768430->r2011": 5,
    "2025768430->r2015": 3,
    "1710476689->r2015": 1,
    "2025768430->r2013": 12,
    "2025768430->r2014": 5,
    "1710476689->r2014": 3,
    "2123284177->r2013": 3,
    "2950320139->r2013": 4,
    "2134842679->r2013": 1,
    "2953267151->r2013": 1,
    "1786904711->r2013": 2,
    "2025768430->r2012": 8,
    "2123284177->r2012": 1,
    "2950320139->r2012": 6,
    "2953267151->r2012": 1,
    "2025768430->r2009": 1,
    "2123284177->r2014": 1,
    "2950320139->r2014": 2,
    "2025768430->r2016": 1,
    "2950320139->r2016": 1,
    "1710476689->r2016": 1,
    "2134842679->r2016": 1,
    "2025768430->r2017": 2,
    "2953267151->r2015": 1,
    "2025768430->r2018": 1,
    "2099471712->r2018": 3,
    "2592298275->r2018": 1,
    "2963865839->r2018": 1,
    "l2004->2123284177": 1,
    "l2007->2025768430": 1,
    "l2012->2099471712": 1,
    "l2009->2025768430": 1,
    "l2009->2099471712": 1,
    "l2009->2123284177": 1,
    "l2009->2159528849": 1,
    "l2009->2421157170": 1,
    "l2009->2953243993": 1,
    "l2009->2963622136": 1,
    "l2013->2099471712": 5,
    "l2013->2134842679": 2,
    "l2013->2159528849": 8,
    "l2013->2953267151": 2,
    "l2013->2953243993": 1,
    "2950320139->r2015": 1,
    "l2014->2592298275": 1,
    "l2014->2751842161": 1,
    "1710476689->r2019": 1,
    "l2014->2134842679": 1,
    "l2013->1843615162": 2,
    "l2005->2123284177": 1,
    "l2005->2134842679": 1,
    "l2005->2953267151": 1,
    "l2013->2770487696": 1,
    "l2006->2123284177": 1,
    "l2006->2134842679": 1,
    "l2006->2953267151": 1,
    "l2015->2421157170": 1,
    "l2015->2963622136": 2,
    "l2015->2159528849": 1,
    "l2015->2593383075": 1,
    "l2014->2159528849": 2,
    "l2014->2770487696": 1,
    "l2014->2953501176": 1,
    "2099471712->r2016": 1,
    "2751842161->r2019": 1,
    "l2016->2593383075": 1,
    "l2016->2953243993": 1,
    "l2016->2980709326": 1,
    "l2016->3093990297": 1,
    "l2018->2953243993": 3,
    "l2019->3001314098": 1,
    "2971040589->r2020": 2,
    "3011483492->r2020": 1,
    "2971040589->r2021": 1,
    "3011483492->r2021": 1
}
const grid = 2;

let leftNodes = [];
let rightNodes = [];

// 解析dot并获取最小年份l和最大年份r
let l = Infinity;
let r = -Infinity;
let labels = '';
let focus_edges_str = '';
let ranks = '';

dot.split('\n').forEach(line => {
    if (line.includes('year')) {
        if (line.includes('rank')) {
            ranks += line + '\n';
        }
    } else if (line.includes('label')) {
        labels += '\t' + line + '\n';
    } else if (line.includes('->')) {
        focus_edges_str += '\t' + line + '\n';
    }
})

// 获取context_edges中 lxxxx 中最小年份 及 rxxxx 中最大年份
Object.keys(context_edges).forEach(edge => {
    let match = edge.match(/l(\d+)->(\d+)/);
    if (match) l = Math.min(l, parseInt(match[1]));
    match = edge.match(/(\d+)->r(\d+)/);
    if (match)  r = Math.max(r, parseInt(match[2]));
});

// 获取ranks中有效的年份（需要除了yearxxxx以外还有其他节点）
let valid_years = []
ranks.split('\n').forEach(line => {
    const match = line.match(/year(\d+) (\d+)/);
    if (match) valid_years.push(parseInt(match[1]));
})

if (l>valid_years[0]) l = valid_years[0];
if (r<valid_years[valid_years.length-1]) r = valid_years[valid_years.length-1];

console.log('valid year range:', l, r)

// 删除ranks中 <l 或 >r 的年份
ranks = ranks.split('\n').filter(line => {
    const match = line.match(/year(\d+)/);
    if(!match) return false;
    const year = parseInt(match[1]);
    return year >= l && year <= r;
}).join('\n');

// 生成每一年的节点连接字符串
for (let year = l; year <= r; year++) {
    leftNodes.push(`l${year}`);
    rightNodes.push(`r${year}`);
}

// 连接生成的节点字符串为一条链
const leftChain = leftNodes.join('->');
const rightChain = rightNodes.join('->');

function replaceYearWithLR(text) {
    // 使用正则表达式匹配 'year' 后跟一串数字，并捕获这些数字
    const regex = /year(\d+)/g;
    // 替换匹配的文本
    return text.replace(regex, (match, year) => `l${year} r${year}`);
}


// 函数用于转换节点名称
function transformNodeName(name, grid) {
    const match = name.match(/^([lr])(\d+)$/);
    if (match) {
        const prefix = match[1];
        const number = parseInt(match[2]);
        if (prefix === 'l') {
            let ret = Math.floor(number / grid) * grid;
            return `l${Math.max(ret, l)}`;
        } else if (prefix === 'r') {
            let ret = (Math.floor(number / grid) + 1) * grid - 1;
            return `r${Math.min(ret, r)}`;
        }
    }
    return name; // 如果没有匹配的lr前缀，返回原始名称
}

// 新的上下文边缘字典
let new_context_edges = {};

// 遍历原始边缘，转换节点名称，并合并边
Object.entries(context_edges).forEach(([edge, weight]) => {
    const parts = edge.split('->');
    const newSrc = transformNodeName(parts[0], grid);
    const newDst = transformNodeName(parts[1], grid);
    const newEdge = `${newSrc}->${newDst}`;

    if (new_context_edges[newEdge]) {
        new_context_edges[newEdge].weight += weight;
        new_context_edges[newEdge].penwidth += weight;  // 这里假设penwidth也是累积的
    } else {
        new_context_edges[newEdge] = { weight: weight, penwidth: weight };
    }
});

// 生成最终的输出字符串
const context_edges_str = Object.entries(new_context_edges).map(([edge, { weight, penwidth }]) => {
    return `${edge} [color="lightgray", weight=${weight}, penwidth=${penwidth}]`;
}).join('\n');

let output = `
digraph G {

crossing_type=0
    
subgraph left {
    style=filled
    color=lightgrey
    node [style=filled,color=lightblue]
    ${leftChain} [weight=10000]
    label = "left"
}

subgraph focus{
    edge [weight=10]
${labels}
${focus_edges_str}
}

subgraph right {
    style=filled
    color=lightgrey
    node [style=filled,color=lightgrey]
    ${rightChain} [weight=10000]
    label = "right"
}

${replaceYearWithLR(ranks)}
${context_edges_str}
}    
`

console.log(output)