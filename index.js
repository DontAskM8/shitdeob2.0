const fs = require("fs")
const beautify = require("js-beautify")
const parser = require("@babel/parser")
const traverse = require("@babel/traverse").default
const generator = require("@babel/generator").default
const types = require("@babel/types")

var original = fs.readFileSync("./raw.js", "utf8")

var undefinedShit = /''\+.*?(?=\))/g
for(let i of original.match(undefinedShit)){
	original = original.replace(i, "undefined")
}


// var evaluableMatches = original.match(/[A-z]+=.*?(,|;)/g)
var evaluableMatches = /(([A-z]|[0-9])+)=(((((\[|\]|\+|!|-| )+)(?=(,|;))))|(((([A-z]|[0-9])+(\*|\+|-|\/)+)+([A-z]|[0-9])+)(?=(,|;))))/g

function blEval(x) {
	var arr = ["{}", "[]", "this", "global"]
	for(let i of arr){
		if((x == i) || (x == (i + ",")) || (x == (i + ";"))){
			return true
		}
	}
	return false
}


// Until line 65 shall be updated to babel but im lazy
var functionsToExecute = /(?<!\w)([A-z]|[0-9])+\(\)/g

for(let k of original.match(functionsToExecute)){
	var contentOfFunction = new RegExp(`(?<=function ${k}\\(\\){).*?(?=})`, "g")
	for(let j of original?.match(contentOfFunction) ?? ""){
		for(let i of j.match(evaluableMatches) ?? ""){
			try {
				var evaled = eval(i)
				if(!blEval(i) && (evaled.toString() != "[object Object]") && !Array.isArray(evaled)){
					eval("global." + i)
					original = original.replace(i, [i.split("=")[0], "=", typeof evaled == "string" ? `'${evaled}'` : evaled].join(""))
				}
			}catch(e){
			}
		}
		
	}
}


for(let i of original.match(evaluableMatches) ?? ""){
	try {
		var evaled = eval(i)
		if(!blEval(i) && (evaled.toString() != "[object Object]") && !Array.isArray(evaled)){
			eval("global." + i)
			original = original.replace(i, [i.split("=")[0], "=", typeof evaled == "string" ? `'${evaled}'` : evaled].join(""))
		}
	}catch(e){
	}
}

var globalVars = Object.keys(global).filter(x => !['global', 'clearInterval', 'clearTimeout', 'setInterval', 'setTimeout', 'queueMicrotask', 'performance', 'clearImmediate', 'setImmediate'].includes(x))
function isGlobalVars(va){
	return globalVars.includes(va) && typeof global[va] != "function"
}

//Main function that triggers the entire universe
var executingFunction = /(?<=(return ))([A-z]|[0-9])+\.call\(this,([A-z]|[0-9])+\)/g

var _testCode = parser.parse(original)
traverse(_testCode, {
	// Replace switch cases var to value
	SwitchStatement: function(path){
		var { node } = path
		node.cases = node.cases.map(x => {
			if(!isGlobalVars(x.test.name)) return x;
			else {
				x.test = types.numericLiteral(global[x.test.name])
				return x
			}
		})
	},
	//Replace var refer var to value
	AssignmentExpression: function(path){
		var { left, right } = path.node
		if(left.type != "Identifier" || right.type != "Identifier" ) return;
		if(!globalVars.includes(right.name)) return;
		try {
			Object.assign(right, types.numericLiteral(global[right.name]))
		}
		catch{}
	},
	//Eval all functions make it global
	FunctionDeclaration: function(path){
		try {
			global[path.node.id.name] = eval(original.slice(path.node.start, path.node.end))
		}catch{
		}
	}
})

original = generator(_testCode, {
	minified: true
}, original).code


var beautifiedCode = beautify(original, {
	indent_size: "\t",
	space_in_empty_paren: true,
	unescape_strings: true
})

fs.writeFileSync("./deobed.js", original)
fs.writeFileSync("./beautified.js", beautifiedCode)